/*===========================================================================
  03_tabelas_resultados.do - Gera tabelas LaTeX da rodada principal
  Andreas Azambuja Barbisan - TCC Insper 2026

  Input:
    checkpoints/estimacao_resumo_tese.csv

  Outputs:
    checkpoints/tabela_principal_eventstudy.tex
    checkpoints/tabela_robustez_sem_imputacao.tex
    checkpoints/tabela_principal_eventstudy.csv
    checkpoints/tabela_robustez_sem_imputacao.csv
===========================================================================*/

clear all
set more off

global root "d:/Andreas/Teste Claude/TCC"
global out  "$root/checkpoints"
global tab  "$root/TCC___Andreas_Azambuja_Barbisan/tabelas"

*===========================================================================
* 1) Carregar resumo da estimacao
*===========================================================================

import delimited "$out/estimacao_resumo_tese.csv", clear varnames(1) stringcols(_all)

* Garantir tipos numericos para campos usados nas tabelas
foreach v in pre_avg pre_sig_count post_avg k0 k1 k2 k3 k4 n_obs n_event_cells {
    capture destring `v', replace force
}

preserve
    use "$out/estimacao_eventstudy_long.dta", clear
    keep if inlist(estimator, "CS21", "DRDID", "SXapr", "TWFE")
    keep if inrange(event_time, 0, 4)
    keep if !missing(att)
    gen double var_event = se^2 if !missing(se)
    collapse (sum) var_event (count) n_post = att, by(serie disc estimator)
    gen double post_se = sqrt(var_event) / n_post if n_post > 0
    keep serie disc estimator post_se
    tempfile post_se
    save `post_se', replace
restore

merge 1:1 serie disc estimator using `post_se', nogen

*===========================================================================
* 2) Tabela principal: media pos-tratamento por estimador
*===========================================================================

preserve
    keep if inlist(estimator, "CS21", "DRDID", "SXapr", "TWFE")
    keep serie disc estimator post_avg post_se pre_avg pre_sig_count n_obs n_event_cells

    reshape wide post_avg post_se pre_avg pre_sig_count n_obs n_event_cells, i(serie disc) j(estimator) string
    sort serie disc

    export delimited using "$out/tabela_principal_eventstudy.csv", replace

    file open ft using "$out/tabela_principal_eventstudy.tex", write replace
    file write ft "\begin{table}[htbp]" _n
    file write ft "\centering" _n
    file write ft "\caption{Efeito médio pós-tratamento (em desvios-padrão)}" _n
    file write ft "\label{tab:eventstudy_principal}" _n
    file write ft "\begin{tabular}{lcccc}" _n
    file write ft "\toprule" _n
    file write ft "Série-disciplina & CS21 & DR-DiD & S\&X aprox. & TWFE \\" _n
    file write ft "\midrule" _n

    forvalues i = 1/`=_N' {
        local cserie = serie[`i']
        local cdisc  = disc[`i']
        local rowlab = "`cserie' `cdisc'"

        foreach e in CS21 DRDID SXapr TWFE {
            if missing(post_avg`e'[`i']) {
                local v_`e' "--"
            }
            else {
                local star ""
                if !missing(post_se`e'[`i']) & post_se`e'[`i'] > 0 {
                    local p = 2 * normal(-abs(post_avg`e'[`i'] / post_se`e'[`i']))
                    if `p' < 0.01 local star "***"
                    else if `p' < 0.05 local star "**"
                    else if `p' < 0.10 local star "*"
                    local v_`e' = "\makecell{" + string(post_avg`e'[`i'], "%6.3f") + "`star'\\(" + string(post_se`e'[`i'], "%6.3f") + ")}"
                }
                else {
                    local v_`e' = string(post_avg`e'[`i'], "%6.3f")
                }
            }
        }

        file write ft "`rowlab' & `v_CS21' & `v_DRDID' & `v_SXapr' & `v_TWFE' \\" _n
    }

    file write ft "\bottomrule" _n
    file write ft "\end{tabular}" _n
    file write ft "\vspace{0.2cm}" _n
    file write ft "\begin{flushleft}" _n
    file write ft "\footnotesize Nota: média de \(k=0\) a \(k=4\) do event study com baseline universal em \(k=-1\). Erros-padrão conservadores entre parênteses, calculados a partir das variâncias dos coeficientes por período; * \(p<0{,}10\), ** \(p<0{,}05\), *** \(p<0{,}01\)." _n
    file write ft "\end{flushleft}" _n
    file write ft "\end{table}" _n
    file close ft
restore

copy "$out/tabela_principal_eventstudy.tex" "$tab/tabela_principal_eventstudy.tex", replace

*===========================================================================
* 3) Tabela robustez: DR-DiD principal vs sem imputacao
*===========================================================================

preserve
    keep if inlist(estimator, "DRDID", "DRRAW")
    keep serie disc estimator post_avg n_obs

    reshape wide post_avg n_obs, i(serie disc) j(estimator) string
    gen n_share = n_obsDRRAW / n_obsDRDID if n_obsDRDID > 0
    sort serie disc

    export delimited using "$out/tabela_robustez_sem_imputacao.csv", replace

    file open fr using "$out/tabela_robustez_sem_imputacao.tex", write replace
    file write fr "\begin{table}[htbp]" _n
    file write fr "\centering" _n
    file write fr "\caption{Sensibilidade à imputação das covariadas}" _n
    file write fr "\label{tab:robustez_sem_imputacao}" _n
    file write fr "\resizebox{\textwidth}{!}{%" _n
    file write fr "\begin{tabular}{lccccc}" _n
    file write fr "\toprule" _n
    file write fr "Série-disciplina & DR-DiD (principal) & DR-DiD sem imputação & N principal & N sem imputação & Fração N \\" _n
    file write fr "\midrule" _n

    forvalues i = 1/`=_N' {
        local cserie = serie[`i']
        local cdisc  = disc[`i']
        local rowlab = "`cserie' `cdisc'"

        local v_main = cond(missing(post_avgDRDID[`i']), "--", string(post_avgDRDID[`i'], "%6.3f"))
        local v_raw  = cond(missing(post_avgDRRAW[`i']), "--", string(post_avgDRRAW[`i'], "%6.3f"))
        local n_main = cond(missing(n_obsDRDID[`i']),    "--", string(n_obsDRDID[`i'],    "%9.0f"))
        local n_raw  = cond(missing(n_obsDRRAW[`i']),    "--", string(n_obsDRRAW[`i'],    "%9.0f"))
        local n_shr  = cond(missing(n_share[`i']),       "--", string(n_share[`i'],       "%6.3f"))

        file write fr "`rowlab' & `v_main' & `v_raw' & `n_main' & `n_raw' & `n_shr' \\" _n
    }

    file write fr "\bottomrule" _n
    file write fr "\end{tabular}" _n
    file write fr "}" _n
    file write fr "\vspace{0.2cm}" _n
    file write fr "\begin{flushleft}" _n
    file write fr "\footnotesize Nota: a versão sem imputação usa apenas observações com covariadas brutas completas. O exercício deve ser lido como diagnóstico de sensibilidade e cobertura amostral, não como substituto direto da especificação principal." _n
    file write fr "\end{flushleft}" _n
    file write fr "\end{table}" _n
    file close fr
restore

copy "$out/tabela_robustez_sem_imputacao.tex" "$tab/tabela_robustez_sem_imputacao.tex", replace

copy "$root/figuras/eventstudy_all_9EF_LP.pdf"  "$tab/imagens/eventstudy_all_9EF_LP.pdf", replace
copy "$root/figuras/eventstudy_all_9EF_Mat.pdf" "$tab/imagens/eventstudy_all_9EF_Mat.pdf", replace
copy "$root/figuras/eventstudy_all_EM3_LP.pdf"  "$tab/imagens/eventstudy_all_EM3_LP.pdf", replace
copy "$root/figuras/eventstudy_all_EM3_Mat.pdf" "$tab/imagens/eventstudy_all_EM3_Mat.pdf", replace

di ""
di "=== Tabelas geradas ==="
di "  - $out/tabela_principal_eventstudy.tex"
di "  - $out/tabela_robustez_sem_imputacao.tex"
di "  - $out/tabela_principal_eventstudy.csv"
di "  - $out/tabela_robustez_sem_imputacao.csv"
di "  - $tab/tabela_principal_eventstudy.tex"
di "  - $tab/tabela_robustez_sem_imputacao.tex"
di "  - $tab/imagens/eventstudy_all_*.pdf"
