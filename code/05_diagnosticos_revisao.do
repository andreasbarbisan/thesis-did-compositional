/*===========================================================================
  05_diagnosticos_revisao.do
  Diagnosticos adicionados apos revisao pre-banca:
    1. Regressao a media: amostra com suporte comparavel de baseline
    2. Imputacao: comparacao entre observacoes completas e imputadas
    3. Spillovers: robustez excluindo controles em regioes com PEI ativo
===========================================================================*/

clear all
set more off

global root "d:/Andreas/Teste Claude/TCC"
global out  "$root/checkpoints"
global tab  "$root/TCC___Andreas_Azambuja_Barbisan/tabelas"

global cov_mun  "ln_pib_pc ln_pop prop_jovem_0_19"
global cov_base "medprof_base $cov_mun"
global cov_comp "pct_fem_imp pct_nee_imp idade_media_imp pct_preta_parda_imp pct_raca_declarada_imp pct_bolsa_fam_imp ln_alunos"
global cov_raw  "pct_fem pct_nee idade_media pct_preta_parda pct_raca_declarada pct_bolsa_fam ln_alunos_raw"

capture mkdir "$out/logs"
capture log close _all
local runstamp = subinstr("`c(current_date)'_`c(current_time)'", " ", "_", .)
local runstamp = subinstr("`runstamp'", ":", "", .)
log using "$out/logs/05_diagnosticos_revisao_`runstamp'.log", text replace

*===========================================================================
* 1. Base de diagnostico
*===========================================================================

use "$out/painel_completo.dta", clear

capture drop medprof_std medprof_year_mean medprof_year_sd
bysort year: egen double medprof_year_mean = mean(medprof)
bysort year: egen double medprof_year_sd   = sd(medprof)
gen double medprof_std = (medprof - medprof_year_mean) / medprof_year_sd ///
    if !missing(medprof) & medprof_year_sd > 0

replace gvar    = 0 if first_treat > 2018 & !missing(first_treat)
replace treated = 0 if first_treat > 2018 & !missing(first_treat)

gen byte cov_completas = !missing(pct_fem, pct_nee, idade_media, ///
    pct_preta_parda, pct_raca_declarada, pct_bolsa_fam, ln_alunos_raw)

capture confirm variable de_id
if _rc encode de, gen(de_id)
capture confirm variable distr_id
if _rc egen distr_id = group(de distr), label

gen byte pei_ativo = first_treat > 0 & first_treat <= year if !missing(first_treat)
replace pei_ativo = 0 if missing(pei_ativo)

bysort year de_id: egen byte de_com_pei_ativo = max(pei_ativo)
bysort year distr_id: egen byte distr_com_pei_ativo = max(pei_ativo)

tempfile base_diag
save `base_diag', replace

*===========================================================================
* 3. Diagnostico de regressao a media
*===========================================================================

tempname preg
postfile `preg' str6 serie str4 disc double drdid_principal drdid_suporte ///
    p10_trat p90_trat n_suporte delta_controles_baixo delta_controles_alto ///
    using "$out/diagnostico_regressao_media.dta", replace

foreach serie in 9EF EM3 {
    foreach disc in 1 2 {
        local dlabel = cond(`disc' == 1, "LP", "Mat")

        use `base_diag', clear
        keep if serie == "`serie'" & co_comp == `disc'
        quietly summarize medprof_base if treated == 1 & !missing(medprof_base), detail
        local p10 = r(p10)
        local p50 = r(p50)
        local p90 = r(p90)

        quietly count if !missing(medprof_std, medprof_base) & ///
            inrange(medprof_base, `p10', `p90')
        local n_suporte = r(N)

        local att_suporte = .
        capture noisily csdid medprof_std $cov_base $cov_comp ///
            if inrange(medprof_base, `p10', `p90'), ///
            ivar(codesc) time(year) gvar(gvar) method(dripw) notyet
        if !_rc {
            capture noisily estat event, window(-4 4) estore(regmedia_`serie'_`disc')
            if !_rc {
                estimates restore regmedia_`serie'_`disc'
                local soma = 0
                local ncoef = 0
                foreach c in Tp0 Tp1 Tp2 Tp3 Tp4 {
                    capture local b = _b[`c']
                    if !_rc & !missing(`b') {
                        local soma = `soma' + `b'
                        local ncoef = `ncoef' + 1
                    }
                }
                if `ncoef' > 0 local att_suporte = `soma' / `ncoef'
            }
        }

        * Convergencia mecanica entre controles nunca tratados, por baseline baixo.
        use `base_diag', clear
        keep if serie == "`serie'" & co_comp == `disc' & gvar == 0
        keep if inlist(year, 2011, 2018) & !missing(medprof_std, medprof_base)
        gen byte baseline_baixo = medprof_base <= `p50' if !missing(medprof_base)
        collapse (mean) medprof_std, by(codesc baseline_baixo year)
        reshape wide medprof_std, i(codesc baseline_baixo) j(year)
        gen delta = medprof_std2018 - medprof_std2011
        quietly summarize delta if baseline_baixo == 1, meanonly
        local delta_low = r(mean)
        quietly summarize delta if baseline_baixo == 0, meanonly
        local delta_high = r(mean)

        local att_main = .
        preserve
            import delimited "$out/estimacao_resumo_tese.csv", clear varnames(1)
            keep if serie == "`serie'" & disc == "`dlabel'" & estimator == "DRDID"
            destring post_avg, replace force
            if _N > 0 local att_main = post_avg[1]
        restore

        post `preg' ("`serie'") ("`dlabel'") (`att_main') (`att_suporte') ///
            (`p10') (`p90') (`n_suporte') (`delta_low') (`delta_high')
    }
}
postclose `preg'

*===========================================================================
* 4. Diagnostico de imputacao
*===========================================================================

use `base_diag', clear
keep if inlist(serie, "9EF", "EM3")

tempname pimp
postfile `pimp' str6 serie str4 disc double n_completa n_imputada ///
    medprof_base_completa medprof_base_imputada tratado_completa tratado_imputada ///
    first_treat_completa first_treat_imputada ///
    using "$out/diagnostico_imputacao.dta", replace

foreach serie in 9EF EM3 {
    foreach disc in 1 2 {
        local dlabel = cond(`disc' == 1, "LP", "Mat")

        preserve
            keep if serie == "`serie'" & co_comp == `disc'
            quietly count if cov_completas == 1
            local n_c = r(N)
            quietly count if cov_completas == 0
            local n_i = r(N)
            quietly summarize medprof_base if cov_completas == 1, meanonly
            local base_c = r(mean)
            quietly summarize medprof_base if cov_completas == 0, meanonly
            local base_i = r(mean)
            quietly summarize treated if cov_completas == 1, meanonly
            local tr_c = r(mean)
            quietly summarize treated if cov_completas == 0, meanonly
            local tr_i = r(mean)
            quietly summarize first_treat if cov_completas == 1 & treated == 1, meanonly
            local ft_c = r(mean)
            quietly summarize first_treat if cov_completas == 0 & treated == 1, meanonly
            local ft_i = r(mean)
        restore

        post `pimp' ("`serie'") ("`dlabel'") (`n_c') (`n_i') ///
            (`base_c') (`base_i') (`tr_c') (`tr_i') (`ft_c') (`ft_i')
    }
}
postclose `pimp'

*===========================================================================
* 5. Robustez de spillover regional
*===========================================================================

tempname pspill
postfile `pspill' str6 serie str4 disc str8 restricao double drdid_principal ///
    drdid_restrito n_restrito using "$out/diagnostico_spillover_regional.dta", replace

foreach serie in 9EF EM3 {
    foreach disc in 1 2 {
        local dlabel = cond(`disc' == 1, "LP", "Mat")

        local att_main = .
        preserve
            import delimited "$out/estimacao_resumo_tese.csv", clear varnames(1)
            keep if serie == "`serie'" & disc == "`dlabel'" & estimator == "DRDID"
            destring post_avg, replace force
            if _N > 0 local att_main = post_avg[1]
        restore

        foreach restr in DE DISTR {
            use `base_diag', clear
            keep if serie == "`serie'" & co_comp == `disc'

            if "`restr'" == "DE" {
                keep if treated == 1 | de_com_pei_ativo == 0
            }
            if "`restr'" == "DISTR" {
                keep if treated == 1 | distr_com_pei_ativo == 0
            }

            quietly count if !missing(medprof_std)
            local n_restr = r(N)
            local att_restr = .

            capture noisily csdid medprof_std $cov_base $cov_comp, ///
                ivar(codesc) time(year) gvar(gvar) method(dripw) notyet
            if !_rc {
                capture noisily estat event, window(-4 4) estore(spill_`restr'_`serie'_`disc')
                if !_rc {
                    estimates restore spill_`restr'_`serie'_`disc'
                    local soma = 0
                    local ncoef = 0
                    foreach c in Tp0 Tp1 Tp2 Tp3 Tp4 {
                        capture local b = _b[`c']
                        if !_rc & !missing(`b') {
                            local soma = `soma' + `b'
                            local ncoef = `ncoef' + 1
                        }
                    }
                    if `ncoef' > 0 local att_restr = `soma' / `ncoef'
                }
            }

            post `pspill' ("`serie'") ("`dlabel'") ("`restr'") ///
                (`att_main') (`att_restr') (`n_restr')
        }
    }
}
postclose `pspill'

*===========================================================================
* 6. Tabelas LaTeX
*===========================================================================

use "$out/diagnostico_regressao_media.dta", clear
export delimited using "$out/diagnostico_regressao_media.csv", replace

file open fr using "$out/tabela_diagnostico_regressao_media.tex", write replace
file write fr "\begin{table}[htbp]" _n
file write fr "\centering" _n
file write fr "\caption{Diagnostico de regressao a media}" _n
file write fr "\label{tab:diagnostico_regressao_media}" _n
file write fr "\resizebox{\textwidth}{!}{%" _n
file write fr "\begin{tabular}{lcccc}" _n
file write fr "\toprule" _n
file write fr "Serie-disciplina & DR-DiD principal & Suporte baseline & $\Delta$ controles baixos & $\Delta$ controles altos \\" _n
file write fr "\midrule" _n
forvalues i = 1/`=_N' {
    local cserie = serie[`i']
    local cdisc = disc[`i']
    local row "`cserie' `cdisc'"
    local a = cond(missing(drdid_principal[`i']), "--", string(drdid_principal[`i'], "%6.3f"))
    local b = cond(missing(drdid_suporte[`i']), "--", string(drdid_suporte[`i'], "%6.3f"))
    local c = cond(missing(delta_controles_baixo[`i']), "--", string(delta_controles_baixo[`i'], "%6.3f"))
    local d = cond(missing(delta_controles_alto[`i']), "--", string(delta_controles_alto[`i'], "%6.3f"))
    file write fr "`row' & `a' & `b' & `c' & `d' \\" _n
}
file write fr "\bottomrule" _n
file write fr "\end{tabular}" _n
file write fr "}" _n
file write fr "\vspace{0.2cm}" _n
file write fr "\begin{flushleft}" _n
file write fr "\footnotesize Nota: a coluna de suporte baseline reestima o DR-DiD restringindo a amostra ao intervalo p10--p90 da proficiencia de linha de base das escolas tratadas. As duas ultimas colunas mostram a variacao media 2011--2018 entre controles nunca tratados abaixo/acima da mediana de baseline das tratadas." _n
file write fr "\end{flushleft}" _n
file write fr "\end{table}" _n
file close fr
copy "$out/tabela_diagnostico_regressao_media.tex" "$tab/tabela_diagnostico_regressao_media.tex", replace

use "$out/diagnostico_imputacao.dta", clear
export delimited using "$out/diagnostico_imputacao.csv", replace

file open fi using "$out/tabela_diagnostico_imputacao.tex", write replace
file write fi "\begin{table}[htbp]" _n
file write fi "\centering" _n
file write fi "\caption{Diagnostico das covariaveis imputadas}" _n
file write fi "\label{tab:diagnostico_imputacao}" _n
file write fi "\begin{tabular}{lcccc}" _n
file write fi "\toprule" _n
file write fi "Serie-disciplina & N completas & N imputadas & Baseline comp. & Baseline imput. \\" _n
file write fi "\midrule" _n
forvalues i = 1/`=_N' {
    local cserie = serie[`i']
    local cdisc = disc[`i']
    local row "`cserie' `cdisc'"
    local nc = string(n_completa[`i'], "%9.0f")
    local ni = string(n_imputada[`i'], "%9.0f")
    local bc = cond(missing(medprof_base_completa[`i']), "--", string(medprof_base_completa[`i'], "%6.2f"))
    local bi = cond(missing(medprof_base_imputada[`i']), "--", string(medprof_base_imputada[`i'], "%6.2f"))
    file write fi "`row' & `nc' & `ni' & `bc' & `bi' \\" _n
}
file write fi "\bottomrule" _n
file write fi "\end{tabular}" _n
file write fi "\vspace{0.2cm}" _n
file write fi "\begin{flushleft}" _n
file write fi "\footnotesize Nota: observacoes completas possuem todas as covariaveis brutas usadas no DR-DiD sem imputacao. Diferencas de baseline entre colunas indicam que a robustez sem imputacao tambem muda a composicao da amostra." _n
file write fi "\end{flushleft}" _n
file write fi "\end{table}" _n
file close fi
copy "$out/tabela_diagnostico_imputacao.tex" "$tab/tabela_diagnostico_imputacao.tex", replace

use "$out/diagnostico_spillover_regional.dta", clear
export delimited using "$out/diagnostico_spillover_regional.csv", replace

file open fs using "$out/tabela_spillover_regional.tex", write replace
file write fs "\begin{table}[htbp]" _n
file write fs "\centering" _n
file write fs "\caption{Robustez a spillovers regionais}" _n
file write fs "\label{tab:spillover_regional}" _n
file write fs "\begin{tabular}{llccc}" _n
file write fs "\toprule" _n
file write fs "Serie-disciplina & Restricao & DR-DiD principal & DR-DiD restrito & N restrito \\" _n
file write fs "\midrule" _n
forvalues i = 1/`=_N' {
    local cserie = serie[`i']
    local cdisc = disc[`i']
    local row "`cserie' `cdisc'"
    local r = restricao[`i']
    local a = cond(missing(drdid_principal[`i']), "--", string(drdid_principal[`i'], "%6.3f"))
    local b = cond(missing(drdid_restrito[`i']), "--", string(drdid_restrito[`i'], "%6.3f"))
    local n = cond(missing(n_restrito[`i']), "--", string(n_restrito[`i'], "%9.0f"))
    file write fs "`row' & `r' & `a' & `b' & `n' \\" _n
}
file write fs "\bottomrule" _n
file write fs "\end{tabular}" _n
file write fs "\vspace{0.2cm}" _n
file write fs "\begin{flushleft}" _n
file write fs "\footnotesize Nota: a restricao remove controles localizados em diretorias de ensino (DE) ou pares DE-distrito (DISTR) que ja tinham ao menos uma escola PEI ativa no ano. E uma proxy regional para contaminacao do contrafactual por spillovers." _n
file write fs "\end{flushleft}" _n
file write fs "\end{table}" _n
file close fs
copy "$out/tabela_spillover_regional.tex" "$tab/tabela_spillover_regional.tex", replace

di ""
di "Diagnosticos de revisao gerados em $out e copiados para $tab."
log close
