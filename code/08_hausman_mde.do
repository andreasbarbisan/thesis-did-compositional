/*===========================================================================
  08_hausman_mde.do
  Regera a tabela de diagnostico tipo Hausman com MDE de 80% de potencia.
===========================================================================*/

clear all
set more off

global root "d:/Andreas/Teste Claude/TCC"
global out  "$root/checkpoints"
global tab  "$root/TCC___Andreas_Azambuja_Barbisan/tabelas"

tempfile resultados drbase

use "$out/estimacao_eventstudy_long.dta", clear
keep if estimator == "DRDID"
keep serie disc event_time att se
rename att att_dr
rename se se_dr
save `drbase', replace

postfile ph str3 serie str3 disc double dr_avg sx_avg diff_avg p_val k_cells mde80 ///
    using `resultados', replace

foreach serie in 9EF EM3 {
foreach disc in LP Mat {
    capture use "$out/sx_`serie'_`disc'.dta", clear
    if _rc continue

    rename l event_time
    keep event_time att se
    rename att att_sx
    rename se se_sx
    gen str3 serie = "`serie'"
    gen str3 disc = "`disc'"

    merge 1:1 serie disc event_time using `drbase', keep(match) nogen
    keep if event_time != -1
    keep if !missing(att_sx, se_sx, att_dr, se_dr) & se_sx > 0 & se_dr > 0

    gen double diff = att_sx - att_dr
    gen double se_diff = sqrt(se_sx^2 + se_dr^2)
    gen double h_piece = (diff / se_diff)^2

    quietly count
    local k = r(N)
    quietly summarize h_piece, meanonly
    local H = r(sum)
    local p = 1 - chi2(`k', `H')

    quietly summarize se_diff, meanonly
    local sebar = r(mean)
    quietly generate double se_diff2 = se_diff^2
    quietly summarize se_diff2, meanonly
    local serms = sqrt(r(mean))
    local mde = (invnormal(0.975) + invnormal(0.80)) * `serms'

    quietly summarize att_dr if inrange(event_time, 0, 4), meanonly
    local dr_avg = r(mean)
    quietly summarize att_sx if inrange(event_time, 0, 4), meanonly
    local sx_avg = r(mean)
    local diff_avg = `sx_avg' - `dr_avg'

    post ph ("`serie'") ("`disc'") (`dr_avg') (`sx_avg') (`diff_avg') (`p') (`k') (`mde')
}
}

postclose ph
use `resultados', clear
export delimited using "$out/tabela_hausman_composicao.csv", replace

file open fh using "$out/tabela_hausman_composicao.tex", write replace
file write fh "\begin{table}[htbp]" _n
file write fh "\centering" _n
file write fh "\caption{Diagnóstico tipo Hausman: DR-DiD versus S\&X aproximado}" _n
file write fh "\label{tab:hausman_composicao}" _n
file write fh "\small" _n
file write fh "\begin{tabular}{lccccc}" _n
file write fh "\toprule" _n
file write fh "Série-disciplina & DR-DiD & S\&X aprox. & Dif. & \(p\)-valor & MDE 80\% \\" _n
file write fh "\midrule" _n

forvalues i = 1/`=_N' {
    local rowlab = cond(serie[`i']=="9EF", "9º Ano EF", "3ª Série EM") + " -- " + cond(disc[`i']=="LP", "LP", "Matemática")
    local dr = subinstr(string(dr_avg[`i'], "%6.3f"), ".", ",", .)
    local sx = subinstr(string(sx_avg[`i'], "%6.3f"), ".", ",", .)
    if diff_avg[`i'] >= 0 {
        local df = "+" + subinstr(string(diff_avg[`i'], "%5.3f"), ".", ",", .)
    }
    else {
        local df = subinstr(string(diff_avg[`i'], "%6.3f"), ".", ",", .)
    }
    local pv = subinstr(string(p_val[`i'], "%6.4f"), ".", ",", .)
    local md = subinstr(string(mde80[`i'], "%6.3f"), ".", ",", .)
    file write fh "`rowlab' & `dr' & `sx' & `df' & `pv' & `md' \\" _n
}

file write fh "\bottomrule" _n
file write fh "\end{tabular}" _n
file write fh "\vspace{0.2cm}" _n
file write fh "\begin{flushleft}" _n
file write fh "\footnotesize Nota: os valores de DR-DiD e S\&X correspondem ao efeito médio pós-tratamento, calculado como a média de \(k=0\) a \(k=4\). O teste usa variância conservadora igual à soma das variâncias marginais. MDE 80\% é a diferença mínima detectável, em desvios-padrão, com 80\% de potência e teste bicaudal de 5\%, calculada a partir da raiz da média dos erros-padrão conservadores ao quadrado nas células comparáveis. A coluna S\&X aprox. exclui células sem suporte comum suficiente e não usa fallback binário para S\&Z." _n
file write fh "\end{flushleft}" _n
file write fh "\end{table}" _n
file close fh

copy "$out/tabela_hausman_composicao.tex" "$tab/tabela_hausman_composicao.tex", replace
