/*===========================================================================
  04_densidade_trat_controle.do - Post-treatment outcome densities
===========================================================================*/

clear all
set more off

global root "d:/Andreas/Teste Claude/TCC"
global out  "$root/checkpoints"
global fig  "$root/figuras"
global tab  "$root/TCC___Andreas_Azambuja_Barbisan/tabelas"
global wpfig "$root/WP_Saresp_Sant_Anna___Xu/files/images"

capture mkdir "$fig"
capture mkdir "$wpfig"
capture mkdir "$out/logs"

capture log close _all
local runstamp = subinstr("`c(current_date)'_`c(current_time)'", " ", "_", .)
local runstamp = subinstr("`runstamp'", ":", "", .)
global density_log "$out/logs/04_densidade_trat_controle_`runstamp'.log"
log using "$density_log", text replace
di "Log desta rodada: $density_log"

use "$out/painel_completo.dta", clear

keep if inlist(serie, "9EF", "EM3") & inlist(co_comp, 1, 2)

capture drop medprof_std medprof_year_mean medprof_year_sd
bysort year: egen double medprof_year_mean = mean(medprof)
bysort year: egen double medprof_year_sd   = sd(medprof)
gen double medprof_std = (medprof - medprof_year_mean) / medprof_year_sd ///
    if !missing(medprof) & medprof_year_sd > 0

gen byte pei_ativo = first_treat > 0 & first_treat <= 2018 & year >= first_treat ///
    if !missing(first_treat)
replace pei_ativo = 0 if missing(pei_ativo)

bysort serie co_comp year: egen byte ano_com_pei_ativo = max(pei_ativo)
gen byte amostra_pos = ano_com_pei_ativo == 1 & !missing(medprof_std)

label define grupo_pos 0 "Control/not-yet" 1 "Active PEI", replace
label values pei_ativo grupo_pos

preserve
    keep if amostra_pos == 1 & !missing(medprof_std)
    collapse (count) N=medprof_std (mean) mean=medprof_std ///
        (sd) sd=medprof_std (p25) p25=medprof_std ///
        (p50) p50=medprof_std (p75) p75=medprof_std, ///
        by(serie co_comp pei_ativo)
    export delimited using "$out/densidade_pos_trat_controle_resumo.csv", replace
restore

local graphs ""

foreach serie in "9EF" "EM3" {
foreach cocomp in 1 2 {

    if `cocomp' == 1 local disc "LP"
    if `cocomp' == 2 local disc "Mat"
    if `cocomp' == 1 local titdisc "Portuguese Language"
    if `cocomp' == 2 local titdisc "Mathematics"
    if "`serie'" == "9EF" local titserie "9th Grade"
    if "`serie'" == "EM3" local titserie "3rd Grade of High School"

    qui count if serie == "`serie'" & co_comp == `cocomp' & ///
        amostra_pos == 1 & pei_ativo == 1 & !missing(medprof_std)
    local n_trat = r(N)
    qui count if serie == "`serie'" & co_comp == `cocomp' & ///
        amostra_pos == 1 & pei_ativo == 0 & !missing(medprof_std)
    local n_ctrl = r(N)

    if `n_trat' > 10 & `n_ctrl' > 10 {
        twoway ///
            (kdensity medprof_std if serie == "`serie'" & co_comp == `cocomp' & ///
                amostra_pos == 1 & pei_ativo == 0, ///
                lcolor(navy) lwidth(medthick)) ///
            (kdensity medprof_std if serie == "`serie'" & co_comp == `cocomp' & ///
                amostra_pos == 1 & pei_ativo == 1, ///
                lcolor(maroon) lwidth(medthick)), ///
            title("`titserie' - `titdisc'", size(medsmall)) ///
            xtitle("Standardized proficiency") ///
            ytitle("Density") ///
            legend(order(1 "Control/not-yet" 2 "Active PEI") ///
                rows(1) ring(0) position(6) size(small) region(style(none))) ///
            graphregion(color(white)) bgcolor(white)

        graph rename dens_`serie'_`disc', replace
        graph export "$fig/densidade_medprof_pos_`serie'_`disc'.pdf", replace
        graph export "$fig/densidade_medprof_pos_`serie'_`disc'.png", width(1800) replace
        copy "$fig/densidade_medprof_pos_`serie'_`disc'.pdf" ///
            "$wpfig/densidade_medprof_pos_`serie'_`disc'.pdf", replace
        local graphs "`graphs' dens_`serie'_`disc'"
    }
    else {
        di "Skipping `serie' `disc': treated N=`n_trat', control N=`n_ctrl'"
    }
}
}

if "`graphs'" != "" {
    graph combine `graphs', cols(2) ///
        title("Density of proficiency in the post-treatment period") ///
        graphregion(color(white))
    graph export "$fig/densidade_medprof_pos_trat_controle.pdf", replace
    graph export "$fig/densidade_medprof_pos_trat_controle.png", width(2400) replace
    copy "$fig/densidade_medprof_pos_trat_controle.pdf" ///
        "$wpfig/densidade_medprof_pos_trat_controle.pdf", replace

    capture copy "$fig/densidade_medprof_pos_trat_controle.pdf" ///
        "$tab/imagens/densidade_medprof_pos_trat_controle.pdf", replace
}

di ""
di "=== Post-treatment density complete ==="
di "Log:      $density_log"
di "Summary:  $out/densidade_pos_trat_controle_resumo.csv"
di "Figures:  $fig/densidade_medprof_pos_trat_controle.pdf"

log close
capture copy "$density_log" "$out/04_densidade_trat_controle_latest.log", replace
