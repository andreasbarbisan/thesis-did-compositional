/*===========================================================================
  06_apendices_estimacoes.do
  Resultados completos para os apÃƒÂªndices empÃƒÂ­ricos

  Objetivo:
    - Organizar outputs por modelo/especificaÃƒÂ§ÃƒÂ£o, nÃƒÂ£o apenas por estimador.
    - Reusar a rodada principal e reestimar diagnÃƒÂ³sticos selecionados.
    - Gerar tabelas completas de ATT por perÃƒÂ­odo e grÃƒÂ¡ficos para o apÃƒÂªndice.
===========================================================================*/

clear all
set more off

global root   "d:/Andreas/Teste Claude/TCC"
global out    "$root/checkpoints"
global fig    "$root/figuras"
global tab    "$root/TCC___Andreas_Azambuja_Barbisan/tabelas"
global apptex "$tab/apendice"
global appfig "$tab/imagens/apendice"
global wpappfig "$root/WP_Saresp_Sant_Anna___Xu/files/images/appendix"

capture mkdir "$apptex"
capture mkdir "$appfig"
capture mkdir "$wpappfig"
capture mkdir "$out/logs"

capture log close _all
local runstamp = subinstr("`c(current_date)'_`c(current_time)'", " ", "_", .)
local runstamp = subinstr("`runstamp'", ":", "", .)
log using "$out/logs/06_apendices_estimacoes_`runstamp'.log", text replace

global K_PRE  4
global K_POST 4

global cov_mun  "ln_pib_pc ln_pop prop_jovem_0_19"
global cov_base "medprof_base $cov_mun"
global cov_comp "pct_fem_imp pct_nee_imp idade_media_imp pct_preta_parda_imp pct_raca_declarada_imp pct_bolsa_fam_imp ln_alunos"

local event_vars ""
forvalues k = $K_PRE(-1)1 {
    local event_vars "`event_vars' Tm`k'"
}
forvalues k = 0/$K_POST {
    local event_vars "`event_vars' Tp`k'"
}

*===========================================================================
* Programas auxiliares
*===========================================================================

capture program drop post_event_rows_model
program define post_event_rows_model
    syntax , Handle(name) Serie(string) Disc(string) Model(string) Eventvars(string)

    capture local nobs = e(N)
    if _rc local nobs = .

    foreach coef of local eventvars {
        if substr("`coef'", 1, 2) == "Tm" {
            local event_time = -real(substr("`coef'", 3, .))
        }
        else {
            local event_time = real(substr("`coef'", 3, .))
        }

        capture local att = _b[`coef']
        if _rc local att = .
        capture local se = _se[`coef']
        if _rc local se = .

        local ci_low = .
        local ci_high = .
        if !missing(`att') & !missing(`se') & `se' > 0 {
            local ci_low  = `att' - 1.96 * `se'
            local ci_high = `att' + 1.96 * `se'
        }

        post `handle' ("`serie'") ("`disc'") ("`model'") ///
            (`event_time') (`att') (`se') (`ci_low') (`ci_high') (`nobs')
    }
end

capture program drop write_period_table_base
program define write_period_table_base
    syntax , Serie(string) Disc(string) File(string) Caption(string) Label(string)

    use "$out/estimacao_eventstudy_long.dta", clear
    keep if serie == "`serie'" & disc == "`disc'"
    keep if inlist(estimator, "TWFE", "CS21", "DRDID", "SXapr")
    keep event_time estimator att se
    reshape wide att se, i(event_time) j(estimator) string
    sort event_time

    file open ft using "`file'", write replace
    file write ft "\begin{table}[!htbp]" _n
    file write ft "\centering" _n
    file write ft "\caption{`caption'}" _n
    file write ft "\label{`label'}" _n
    file write ft "\vspace{-0.15cm}" _n
    file write ft "\scriptsize" _n
    file write ft "\setlength{\tabcolsep}{4pt}" _n
    file write ft "\renewcommand{\arraystretch}{0.92}" _n
    file write ft "\resizebox{\textwidth}{!}{%" _n
    file write ft "\begin{tabular}{lcccc}" _n
    file write ft "\toprule" _n
    file write ft "PerÃƒÂ­odo relativo & TWFE & CS21 & DR-DiD & S\&X aprox. \\" _n
    file write ft "\midrule" _n

    forvalues i = 1/`=_N' {
        local kk = string(event_time[`i'], "%9.0f")
        foreach e in TWFE CS21 DRDID SXapr {
            if missing(att`e'[`i']) {
                local c_`e' "--"
            }
            else if missing(se`e'[`i']) | se`e'[`i'] == 0 {
                local c_`e' = string(att`e'[`i'], "%6.3f")
            }
            else {
                local c_`e' = string(att`e'[`i'], "%6.3f") + " (" + string(se`e'[`i'], "%6.3f") + ")"
            }
        }
        file write ft "`kk' & `c_TWFE' & `c_CS21' & `c_DRDID' & `c_SXapr' \\" _n
    }

    file write ft "\bottomrule" _n
    file write ft "\end{tabular}" _n
    file write ft "}" _n
    file write ft "\begin{flushleft}" _n
    file write ft "\footnotesize Nota: cada cÃƒÂ©lula reporta ATT e erro-padrÃƒÂ£o entre parÃƒÂªnteses, em desvios-padrÃƒÂ£o da proficiÃƒÂªncia anual. O perÃƒÂ­odo \(k=-1\) ÃƒÂ© normalizado como linha de base; por isso aparece sem erro-padrÃƒÂ£o." _n
    file write ft "\end{flushleft}" _n
    file write ft "\end{table}" _n
    file close ft
end

capture program drop write_period_table_models
program define write_period_table_models
    syntax , Serie(string) Disc(string) File(string) Caption(string) Label(string) Models(string) Header(string)

    use "$out/apendice_model_eventstudy_long.dta", clear
    keep if serie == "`serie'" & disc == "`disc'"
    gen byte keep_model = 0
    foreach m of local models {
        replace keep_model = 1 if model == "`m'"
    }
    keep if keep_model == 1
    keep event_time model att se
    reshape wide att se, i(event_time) j(model) string
    sort event_time

    local nmodels : word count `models'
    local colspec "l"
    forvalues j = 1/`nmodels' {
        local colspec "`colspec'c"
    }

    file open ft using "`file'", write replace
    file write ft "\begin{table}[!htbp]" _n
    file write ft "\centering" _n
    file write ft "\caption{`caption'}" _n
    file write ft "\label{`label'}" _n
    file write ft "\vspace{-0.15cm}" _n
    file write ft "\scriptsize" _n
    file write ft "\setlength{\tabcolsep}{4pt}" _n
    file write ft "\renewcommand{\arraystretch}{0.92}" _n
    file write ft "\begin{tabular}{`colspec'}" _n
    file write ft "\toprule" _n
    file write ft "PerÃƒÂ­odo relativo & `header' \\" _n
    file write ft "\midrule" _n

    forvalues i = 1/`=_N' {
        local kk = string(event_time[`i'], "%9.0f")
        local row "`kk'"
        foreach m of local models {
            capture confirm variable att`m'
            if _rc {
                local cell "--"
            }
            else if missing(att`m'[`i']) {
                local cell "--"
            }
            else if missing(se`m'[`i']) | se`m'[`i'] == 0 {
                local cell = string(att`m'[`i'], "%6.3f")
            }
            else {
                local cell = string(att`m'[`i'], "%6.3f") + " (" + string(se`m'[`i'], "%6.3f") + ")"
            }
            local row "`row' & `cell'"
        }
        file write ft "`row' \\" _n
    }

    file write ft "\bottomrule" _n
    file write ft "\end{tabular}" _n
    file write ft "\begin{flushleft}" _n
    file write ft "\footnotesize Nota: cada cÃƒÂ©lula reporta ATT e erro-padrÃƒÂ£o entre parÃƒÂªnteses, em desvios-padrÃƒÂ£o da proficiÃƒÂªncia anual. O perÃƒÂ­odo \(k=-1\) ÃƒÂ© normalizado como linha de base; por isso aparece sem erro-padrÃƒÂ£o." _n
    file write ft "\end{flushleft}" _n
    file write ft "\end{table}" _n
    file close ft
end

*===========================================================================
* 1. Tabelas completas do modelo principal
*===========================================================================

copy "$fig/eventstudy_all_9EF_LP.pdf"  "$appfig/eventstudy_base_9EF_LP.pdf", replace
copy "$fig/eventstudy_all_9EF_LP.pdf"  "$wpappfig/eventstudy_base_9EF_LP.pdf", replace
copy "$fig/eventstudy_all_9EF_Mat.pdf" "$appfig/eventstudy_base_9EF_Mat.pdf", replace
copy "$fig/eventstudy_all_9EF_Mat.pdf" "$wpappfig/eventstudy_base_9EF_Mat.pdf", replace
copy "$fig/eventstudy_all_EM3_LP.pdf"  "$appfig/eventstudy_base_EM3_LP.pdf", replace
copy "$fig/eventstudy_all_EM3_LP.pdf"  "$wpappfig/eventstudy_base_EM3_LP.pdf", replace
copy "$fig/eventstudy_all_EM3_Mat.pdf" "$appfig/eventstudy_base_EM3_Mat.pdf", replace
copy "$fig/eventstudy_all_EM3_Mat.pdf" "$wpappfig/eventstudy_base_EM3_Mat.pdf", replace

write_period_table_base, serie("9EF") disc("LP") ///
    file("$apptex/ap_event_base_9EF_LP.tex") ///
    caption("Modelo principal, 9Ã‚Âº ano EF em LÃƒÂ­ngua Portuguesa: ATT por perÃƒÂ­odo relativo") ///
    label("tab:ap_event_base_9ef_lp")

write_period_table_base, serie("9EF") disc("Mat") ///
    file("$apptex/ap_event_base_9EF_Mat.tex") ///
    caption("Modelo principal, 9Ã‚Âº ano EF em MatemÃƒÂ¡tica: ATT por perÃƒÂ­odo relativo") ///
    label("tab:ap_event_base_9ef_mat")

write_period_table_base, serie("EM3") disc("LP") ///
    file("$apptex/ap_event_base_EM3_LP.tex") ///
    caption("Modelo principal, 3Ã‚Âª sÃƒÂ©rie EM em LÃƒÂ­ngua Portuguesa: ATT por perÃƒÂ­odo relativo") ///
    label("tab:ap_event_base_em3_lp")

write_period_table_base, serie("EM3") disc("Mat") ///
    file("$apptex/ap_event_base_EM3_Mat.tex") ///
    caption("Modelo principal, 3Ã‚Âª sÃƒÂ©rie EM em MatemÃƒÂ¡tica: ATT por perÃƒÂ­odo relativo") ///
    label("tab:ap_event_base_em3_mat")

*===========================================================================
* 2. Reestimacoes por modelo: regressao a media e spillovers regionais
*===========================================================================

use "$out/estimacao_eventstudy_long.dta", clear
keep if estimator == "DRDID" | estimator == "DRRAW"
gen str12 model = cond(estimator == "DRDID", "Principal", "SemImput")
keep serie disc model event_time att se ci_low ci_high N_obs
rename N_obs n_obs
tempfile modelos_existentes
save `modelos_existentes', replace

use "$out/painel_completo.dta", clear

capture drop medprof_std medprof_year_mean medprof_year_sd
bysort year: egen double medprof_year_mean = mean(medprof)
bysort year: egen double medprof_year_sd   = sd(medprof)
gen double medprof_std = (medprof - medprof_year_mean) / medprof_year_sd ///
    if !missing(medprof) & medprof_year_sd > 0

replace gvar    = 0 if first_treat > 2018 & !missing(first_treat)
replace treated = 0 if first_treat > 2018 & !missing(first_treat)

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

tempname pevent
postfile `pevent' str6 serie str4 disc str12 model int event_time ///
    double att se ci_low ci_high n_obs using "$out/apendice_model_eventstudy_extra.dta", replace

foreach serie in 9EF EM3 {
    foreach discnum in 1 2 {
        local dlabel = cond(`discnum' == 1, "LP", "Mat")

        use `base_diag', clear
        keep if serie == "`serie'" & co_comp == `discnum'
        quietly summarize medprof_base if treated == 1 & !missing(medprof_base), detail
        local p10 = r(p10)
        local p90 = r(p90)

        capture noisily csdid medprof_std $cov_base $cov_comp ///
            if inrange(medprof_base, `p10', `p90'), ///
            ivar(codesc) time(year) gvar(gvar) method(dripw) notyet
        if !_rc {
            capture noisily estat event, window(-4 4) estore(ap_RegMedia_`serie'_`discnum')
            if !_rc {
                estimates restore ap_RegMedia_`serie'_`discnum'
                post_event_rows_model, handle(`pevent') serie("`serie'") ///
                    disc("`dlabel'") model("RegMedia") eventvars("`event_vars'")
            }
        }

        foreach restr in DE DISTR {
            use `base_diag', clear
            keep if serie == "`serie'" & co_comp == `discnum'

            if "`restr'" == "DE" {
                keep if treated == 1 | de_com_pei_ativo == 0
                local modelname "SpillDE"
            }
            if "`restr'" == "DISTR" {
                keep if treated == 1 | distr_com_pei_ativo == 0
                local modelname "SpillDISTR"
            }

            capture noisily csdid medprof_std $cov_base $cov_comp, ///
                ivar(codesc) time(year) gvar(gvar) method(dripw) notyet
            if !_rc {
                capture noisily estat event, window(-4 4) estore(ap_`modelname'_`serie'_`discnum')
                if !_rc {
                    estimates restore ap_`modelname'_`serie'_`discnum'
                    post_event_rows_model, handle(`pevent') serie("`serie'") ///
                        disc("`dlabel'") model("`modelname'") eventvars("`event_vars'")
                }
            }
        }
    }
}

postclose `pevent'

use `modelos_existentes', clear
append using "$out/apendice_model_eventstudy_extra.dta"
sort serie disc model event_time
save "$out/apendice_model_eventstudy_long.dta", replace
export delimited using "$out/apendice_model_eventstudy_long.csv", replace

*===========================================================================
* 3. Graficos comparando modelos
*===========================================================================

foreach serie in 9EF EM3 {
    foreach disc in LP Mat {
        local titserie = cond("`serie'" == "9EF", "9th Grade", "3rd Grade of High School")
        local titdisc  = cond("`disc'" == "LP", "Portuguese Language", "Mathematics")

        use "$out/apendice_model_eventstudy_long.dta", clear
        keep if serie == "`serie'" & disc == "`disc'" & inlist(model, "Principal", "RegMedia")
        twoway ///
            (rcap ci_low ci_high event_time if model == "Principal", lcolor(navy%60)) ///
            (connected att event_time if model == "Principal", sort lcolor(navy) mcolor(navy) msymbol(O)) ///
            (rcap ci_low ci_high event_time if model == "RegMedia", lcolor(maroon%60)) ///
            (connected att event_time if model == "RegMedia", sort lcolor(maroon) mcolor(maroon) msymbol(D)), ///
            yline(0, lcolor(black) lpattern(dash)) ///
            xline(-0.5, lcolor(gs8) lpattern(dash)) ///
            xlabel(-4(1)4) ///
            title("`titserie' - `titdisc'") ///
            subtitle("Main specification vs. P10-P90 support") ///
            xtitle("Relative time") ///
            ytitle("ATT") ///
            legend(order(2 "Main specification" 4 "P10-P90 support") rows(1) ring(0) position(6) size(small) region(style(none))) ///
            graphregion(color(white)) bgcolor(white)
        graph export "$appfig/eventstudy_regmedia_`serie'_`disc'.pdf", replace
        copy "$appfig/eventstudy_regmedia_`serie'_`disc'.pdf" "$wpappfig/eventstudy_regmedia_`serie'_`disc'.pdf", replace

        use "$out/apendice_model_eventstudy_long.dta", clear
        keep if serie == "`serie'" & disc == "`disc'" & inlist(model, "Principal", "SpillDE", "SpillDISTR")
        twoway ///
            (rcap ci_low ci_high event_time if model == "Principal", lcolor(navy%55)) ///
            (connected att event_time if model == "Principal", sort lcolor(navy) mcolor(navy) msymbol(O)) ///
            (rcap ci_low ci_high event_time if model == "SpillDE", lcolor(forest_green%55)) ///
            (connected att event_time if model == "SpillDE", sort lcolor(forest_green) mcolor(forest_green) msymbol(T)) ///
            (rcap ci_low ci_high event_time if model == "SpillDISTR", lcolor(maroon%55)) ///
            (connected att event_time if model == "SpillDISTR", sort lcolor(maroon) mcolor(maroon) msymbol(D)), ///
            yline(0, lcolor(black) lpattern(dash)) ///
            xline(-0.5, lcolor(gs8) lpattern(dash)) ///
            xlabel(-4(1)4) ///
            title("`titserie' - `titdisc'") ///
            subtitle("Main specification vs. regional filter") ///
            xtitle("Relative time") ///
            ytitle("ATT") ///
            legend(order(2 "Main specification" 4 "Exclude exposed DE" 6 "Exclude exposed district-DE") ///
                rows(1) ring(0) position(6) size(small) region(style(none))) ///
            graphregion(color(white)) bgcolor(white)
        graph export "$appfig/eventstudy_spillover_`serie'_`disc'.pdf", replace
        copy "$appfig/eventstudy_spillover_`serie'_`disc'.pdf" "$wpappfig/eventstudy_spillover_`serie'_`disc'.pdf", replace

        use "$out/apendice_model_eventstudy_long.dta", clear
        keep if serie == "`serie'" & disc == "`disc'" & inlist(model, "Principal", "SemImput")
        quietly count if model == "SemImput" & !missing(att)
        if r(N) > 0 {
            twoway ///
                (rcap ci_low ci_high event_time if model == "Principal", lcolor(navy%60)) ///
                (connected att event_time if model == "Principal", sort lcolor(navy) mcolor(navy) msymbol(O)) ///
                (rcap ci_low ci_high event_time if model == "SemImput", lcolor(maroon%60)) ///
                (connected att event_time if model == "SemImput", sort lcolor(maroon) mcolor(maroon) msymbol(D)), ///
                yline(0, lcolor(black) lpattern(dash)) ///
                xline(-0.5, lcolor(gs8) lpattern(dash)) ///
                xlabel(-4(1)4) ///
                title("`titserie' - `titdisc'") ///
                subtitle("Main specification vs. complete-case sample") ///
                xtitle("Relative time") ///
                ytitle("ATT") ///
                legend(order(2 "Main specification" 4 "No imputation") rows(1) ring(0) position(6) size(small) region(style(none))) ///
                graphregion(color(white)) bgcolor(white)
            graph export "$appfig/eventstudy_imputacao_`serie'_`disc'.pdf", replace
            copy "$appfig/eventstudy_imputacao_`serie'_`disc'.pdf" "$wpappfig/eventstudy_imputacao_`serie'_`disc'.pdf", replace
        }
    }
}

di "Appendix figures regenerated. Tables were intentionally left untouched."
log close
exit, clear

*===========================================================================
* 4. Tabelas completas por modelo de diagnostico
*===========================================================================

foreach serie in 9EF EM3 {
    foreach disc in LP Mat {
        local lbls = lower("`serie'_`disc'")
        local lbls = subinstr("`lbls'", "mat", "mat", .)

        write_period_table_models, serie("`serie'") disc("`disc'") ///
            file("$apptex/ap_event_regmedia_`serie'_`disc'.tex") ///
            caption("DiagnÃƒÂ³stico de regressÃƒÂ£o ÃƒÂ  mÃƒÂ©dia, `serie' `disc': ATT por perÃƒÂ­odo relativo") ///
            label("tab:ap_event_regmedia_`lbls'") ///
            models("Principal RegMedia") ///
            header("Principal & Suporte P10-P90")

        write_period_table_models, serie("`serie'") disc("`disc'") ///
            file("$apptex/ap_event_spillover_`serie'_`disc'.tex") ///
            caption("Robustez a spillovers regionais, `serie' `disc': ATT por perÃƒÂ­odo relativo") ///
            label("tab:ap_event_spillover_`lbls'") ///
            models("Principal SpillDE SpillDISTR") ///
            header("Principal & Exclui DE exposta & Exclui distrito-DE exposto")

        write_period_table_models, serie("`serie'") disc("`disc'") ///
            file("$apptex/ap_event_imputacao_`serie'_`disc'.tex") ///
            caption("DiagnÃƒÂ³stico de imputaÃƒÂ§ÃƒÂ£o, `serie' `disc': ATT por perÃƒÂ­odo relativo") ///
            label("tab:ap_event_imputacao_`lbls'") ///
            models("Principal SemImput") ///
            header("Principal & Sem imputaÃƒÂ§ÃƒÂ£o")
    }
}

*===========================================================================
* 5. Detalhes S&X: celulas por coorte e periodo relativo
*===========================================================================

use "$out/estimacao_sx_detalhe.dta", clear
gen double variance = se^2 if !missing(se)
export delimited using "$out/apendice_sx_celulas.csv", replace

keep serie disc rel_time cohort cal_time n_sample n_treated att se variance status
sort serie disc rel_time cohort

foreach serie in 9EF EM3 {
    foreach disc in LP Mat {
        preserve
            keep if serie == "`serie'" & disc == "`disc'"
                local slower = lower("`serie'")
                local dlower = lower("`disc'")

                file open fsx using "$apptex/ap_sx_celulas_`serie'_`disc'.tex", write replace
                file write fsx "\begingroup" _n
                file write fsx "\tiny" _n
                file write fsx "\setlength{\tabcolsep}{2.2pt}" _n
                file write fsx "\renewcommand{\arraystretch}{0.82}" _n
                file write fsx "\begin{longtable}{rrrrrrrll}" _n
                file write fsx "\caption{AproximaÃƒÂ§ÃƒÂ£o S\&X, `serie' `disc': cÃƒÂ©lulas coorte-perÃƒÂ­odo usadas na agregaÃƒÂ§ÃƒÂ£o}" _n
                file write fsx "\label{tab:ap_sx_celulas_`slower'_`dlower'}\\" _n
                file write fsx "\toprule" _n
                file write fsx "\(k\) & Coorte & Ano & N & N trat. & ATT & EP & Var. & Uso \\" _n
                file write fsx "\midrule" _n
                file write fsx "\endfirsthead" _n
                file write fsx "\toprule" _n
                file write fsx "\(k\) & Coorte & Ano & N & N trat. & ATT & EP & Var. & Uso \\" _n
                file write fsx "\midrule" _n
                file write fsx "\endhead" _n
                file write fsx "\midrule" _n
                file write fsx "\multicolumn{9}{r}{Continua na prÃƒÂ³xima pÃƒÂ¡gina} \\" _n
                file write fsx "\endfoot" _n
                file write fsx "\bottomrule" _n
                file write fsx "\endlastfoot" _n

                forvalues i = 1/`=_N' {
                    local k  = string(rel_time[`i'], "%4.0f")
                    local g  = string(cohort[`i'], "%4.0f")
                    local yy = string(cal_time[`i'], "%4.0f")
                    local nn = string(n_sample[`i'], "%9.0f")
                    local nt = string(n_treated[`i'], "%9.0f")
                    local a  = cond(missing(att[`i']), "--", string(att[`i'], "%6.3f"))
                    local s  = cond(missing(se[`i']), "--", string(se[`i'], "%6.3f"))
                    local v  = cond(missing(variance[`i']), "--", string(variance[`i'], "%6.4f"))
                    local st_raw = status[`i']
                    local st "Estimado"
                    if "`st_raw'" == "ps_lowoverlap_k5" local st "Sem suporte"
                    if "`st_raw'" == "cell_too_small" local st "CÃƒÂ©lula pequena"
                    if "`st_raw'" == "baseline_km1" local st "Base"
                    file write fsx "`k' & `g' & `yy' & `nn' & `nt' & `a' & `s' & `v' & `st' \\" _n
                }

                file write fsx "\end{longtable}" _n
                file write fsx "\vspace{-0.35cm}" _n
                file write fsx "\noindent\footnotesize Nota: \textit{Estimado} indica cÃƒÂ©lula usada na mÃƒÂ©dia S\&X; \textit{Sem suporte} indica baixa sobreposiÃƒÂ§ÃƒÂ£o no logit multinomial; \textit{CÃƒÂ©lula pequena} indica ausÃƒÂªncia ou nÃƒÂºmero insuficiente de tratados; \textit{Base} ÃƒÂ© o perÃƒÂ­odo normalizado \(k=-1\)." _n
                file write fsx "\par\endgroup" _n
                file close fsx
        restore
    }
}

preserve
    collapse (count) n_celulas=att (mean) att_medio=att se_medio=se variance_media=variance, by(serie disc rel_time)
    sort serie disc rel_time

    file open fsum using "$apptex/ap_sx_resumo_periodo.tex", write replace
    file write fsum "\begin{table}[!htbp]" _n
    file write fsum "\centering" _n
    file write fsum "\caption{AproximaÃƒÂ§ÃƒÂ£o S\&X: resumo das cÃƒÂ©lulas estimadas por perÃƒÂ­odo relativo}" _n
    file write fsum "\label{tab:ap_sx_resumo_periodo}" _n
    file write fsum "\vspace{-0.15cm}" _n
    file write fsum "\scriptsize" _n
    file write fsum "\setlength{\tabcolsep}{4pt}" _n
    file write fsum "\renewcommand{\arraystretch}{0.92}" _n
    file write fsum "\begin{tabular}{llrrrr}" _n
    file write fsum "\toprule" _n
    file write fsum "SÃƒÂ©rie-disciplina & \(k\) & CÃƒÂ©lulas & ATT mÃƒÂ©dio & EP mÃƒÂ©dio & Var. mÃƒÂ©dia \\" _n
    file write fsum "\midrule" _n

    forvalues i = 1/`=_N' {
        local row = serie[`i'] + " " + disc[`i']
        local k = string(rel_time[`i'], "%4.0f")
        local n = string(n_celulas[`i'], "%4.0f")
        local a = cond(missing(att_medio[`i']), "--", string(att_medio[`i'], "%6.3f"))
        local s = cond(missing(se_medio[`i']), "--", string(se_medio[`i'], "%6.3f"))
        local v = cond(missing(variance_media[`i']), "--", string(variance_media[`i'], "%6.4f"))
        file write fsum "`row' & `k' & `n' & `a' & `s' & `v' \\" _n
    }

    file write fsum "\bottomrule" _n
    file write fsum "\end{tabular}" _n
    file write fsum "\begin{flushleft}" _n
    file write fsum "\footnotesize Nota: o nÃƒÂºmero de cÃƒÂ©lulas corresponde ÃƒÂ s combinaÃƒÂ§ÃƒÂµes coorte-perÃƒÂ­odo com suporte suficiente para entrar na mÃƒÂ©dia S\&X aproximada. ATT, erro-padrÃƒÂ£o e variÃƒÂ¢ncia sÃƒÂ£o mÃƒÂ©dias simples entre as cÃƒÂ©lulas estimadas naquele perÃƒÂ­odo relativo." _n
    file write fsum "\end{flushleft}" _n
    file write fsum "\end{table}" _n
    file close fsum
restore

di _newline "=== Apendices empiricos gerados ==="
di "Tabelas: $apptex"
di "Figuras: $appfig"
di "Dados:   $out/apendice_model_eventstudy_long.csv"

log close
