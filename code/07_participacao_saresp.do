/*===========================================================================
  07_participacao_saresp.do
  Diagnostico simples da taxa de participacao no SARESP

  Numerador: alunos com participacao registrada na prova, por escola-ano-serie.
  Denominador: matriculas estaduais ativas por escola-ano-serie.
===========================================================================*/

clear all
set more off

global root  "d:/Andreas/Teste Claude/TCC"
global micro "$root/Microdados Saresp"
global out   "$root/checkpoints"
global tab   "$root/TCC___Andreas_Azambuja_Barbisan/tabelas"

capture mkdir "$out/logs"
capture log close _all
local runstamp = subinstr("`c(current_date)'_`c(current_time)'", " ", "_", .)
local runstamp = subinstr("`runstamp'", ":", "", .)
log using "$out/logs/07_participacao_saresp_`runstamp'.log", text replace

local f2011 "$micro/MICRODADOS_SARESP_2011_0.csv"
local f2012 "$micro/MICRODADOS_SARESP_2012.csv"
local f2013 "$micro/MICRODADOS_SARESP_2013.csv"
local f2014 "$micro/MICRODADOS_SARESP_2014.csv"
local f2015 "$micro/MICRODADOS SARESP 2015 - DADOS ABERTO.csv"
local f2016 "$micro/MICRODADOS SARESP 2016 - DADOS ABERTO.csv"
local f2017 "$micro/MICRODADOS SARESP 2017 - DADOS ABERTO.csv"
local f2018 "$micro/MICRODADOS SARESP 2018 - DADOS ABERTO_0.csv"

tempfile base_part
local first = 1

foreach ano in 2011 2012 2013 2014 2015 2016 2017 2018 {
    import delimited "`f`ano''", delimiter(";") encoding("latin1") ///
        clear case(lower) varnames(1)

    capture rename serie_ano serie_raw
    capture rename serie serie_raw
    tostring serie_raw, replace force
    gen str3 serie = ""
    replace serie = "9EF" if strpos(serie_raw, "9") > 0
    replace serie = "EM3" if serie_raw == "3" | strpos(upper(serie_raw), "EM-3") > 0 | ///
        (strpos(serie_raw, "3") > 0 & strpos(upper(serie_raw), "EM") > 0)
    keep if inlist(serie, "9EF", "EM3")

    destring codesc particip_lp particip_mat, replace force
    gen year = `ano'

    preserve
        gen str3 disc = "LP"
        gen byte participou = particip_lp == 1 if !missing(particip_lp)
        replace participou = 0 if missing(participou)
        collapse (sum) n_particip = participou (count) n_reg_saresp = participou, ///
            by(codesc year serie disc)
        tempfile lp`ano'
        save `lp`ano'', replace
    restore

    preserve
        gen str3 disc = "Mat"
        gen byte participou = particip_mat == 1 if !missing(particip_mat)
        replace participou = 0 if missing(participou)
        collapse (sum) n_particip = participou (count) n_reg_saresp = participou, ///
            by(codesc year serie disc)
        append using `lp`ano''
        tempfile ano`ano'
        save `ano`ano'', replace
    restore

    if `first' {
        use `ano`ano'', clear
        save `base_part', replace
        local first = 0
    }
    else {
        use `base_part', clear
        append using `ano`ano''
        save `base_part', replace
    }
}

use `base_part', clear

preserve
    use "$out/painel_completo.dta", clear
    keep codesc year serie n_alunos_imp first_treat
    duplicates drop
    rename n_alunos_imp n_alunos
    tempfile denom
    save `denom', replace
restore

merge m:1 codesc year serie using `denom', keep(master match) ///
    keepusing(n_alunos first_treat) nogen
replace first_treat = 0 if missing(first_treat)

gen double taxa_participacao = n_particip / n_alunos if n_alunos > 0
replace taxa_participacao = 1 if taxa_participacao > 1 & !missing(taxa_participacao)

gen byte status = .
replace status = 0 if first_treat == 0
replace status = 1 if first_treat > 0 & year < first_treat
replace status = 2 if first_treat > 0 & year >= first_treat
label define status_part 0 "Controle" 1 "Tratado pre" 2 "Tratado pos"
label values status status_part

save "$out/diagnostico_participacao_saresp.dta", replace
export delimited using "$out/diagnostico_participacao_saresp.csv", replace

preserve
    keep if !missing(taxa_participacao) & inlist(status, 0, 1, 2)
    collapse (mean) taxa_participacao (sum) n_particip n_alunos ///
        (count) n_escola_ano = taxa_participacao, by(serie disc status)
    reshape wide taxa_participacao n_particip n_alunos n_escola_ano, ///
        i(serie disc) j(status)

    gen double diff_pos_pre = taxa_participacao2 - taxa_participacao1
    gen double n_share_pos_pre = n_escola_ano2 / n_escola_ano1 if n_escola_ano1 > 0
    sort serie disc

    export delimited using "$out/tabela_participacao_saresp.csv", replace

    file open fp using "$out/tabela_participacao_saresp.tex", write replace
    file write fp "\begin{table}[htbp]" _n
    file write fp "\centering" _n
    file write fp "\caption{Diagnóstico da taxa de participação no SARESP}" _n
    file write fp "\label{tab:participacao_saresp}" _n
    file write fp "\small" _n
    file write fp "\begin{tabular}{lcccc}" _n
    file write fp "\toprule" _n
    file write fp "Série-disciplina & Controle & Tratado pré & Tratado pós & Dif. pós--pré \\" _n
    file write fp "\midrule" _n

    forvalues i = 1/`=_N' {
        local rowlab = serie[`i'] + " " + disc[`i']
        local c0 = cond(missing(taxa_participacao0[`i']), "--", string(100 * taxa_participacao0[`i'], "%6.1f") + "\%")
        local c1 = cond(missing(taxa_participacao1[`i']), "--", string(100 * taxa_participacao1[`i'], "%6.1f") + "\%")
        local c2 = cond(missing(taxa_participacao2[`i']), "--", string(100 * taxa_participacao2[`i'], "%6.1f") + "\%")
        if missing(diff_pos_pre[`i']) {
            local cd "--"
        }
        else if diff_pos_pre[`i'] >= 0 {
            local cd = "+" + string(100 * diff_pos_pre[`i'], "%4.1f") + " p.p."
        }
        else {
            local cd = string(100 * diff_pos_pre[`i'], "%5.1f") + " p.p."
        }
        file write fp "`rowlab' & `c0' & `c1' & `c2' & `cd' \\" _n
    }

    file write fp "\bottomrule" _n
    file write fp "\end{tabular}" _n
    file write fp "\vspace{0.2cm}" _n
    file write fp "\begin{flushleft}" _n
    file write fp "\footnotesize Nota: a taxa divide o número de alunos com participação registrada na prova pelo total de matrículas estaduais ativas na escola-série-ano, usando o denominador harmonizado do painel final. A diferença pós--pré é descritiva e restrita às escolas que aderem ao PEI entre 2012 e 2018." _n
    file write fp "\end{flushleft}" _n
    file write fp "\end{table}" _n
    file close fp
restore

copy "$out/tabela_participacao_saresp.tex" "$tab/tabela_participacao_saresp.tex", replace

di "Tabela salva em $tab/tabela_participacao_saresp.tex"
log close
