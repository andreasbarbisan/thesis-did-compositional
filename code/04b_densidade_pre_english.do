/*===========================================================================
  04b_densidade_pre_english.do - Pre-treatment outcome density figure
===========================================================================*/

clear all
set more off

global root "d:/Andreas/Teste Claude/TCC"
global out  "$root/checkpoints"
global fig  "$root/figuras"
global wpfig "$root/WP_Saresp_Sant_Anna___Xu/files/images"

capture mkdir "$fig"
capture mkdir "$wpfig"

use "$out/painel_completo.dta", clear
keep if inlist(serie, "9EF", "EM3") & inlist(co_comp, 1, 2)

capture drop medprof_std medprof_year_mean medprof_year_sd
bysort year: egen double medprof_year_mean = mean(medprof)
bysort year: egen double medprof_year_sd   = sd(medprof)
gen double medprof_std = (medprof - medprof_year_mean) / medprof_year_sd ///
    if !missing(medprof) & medprof_year_sd > 0

gen byte pre_treated = first_treat > 0 & first_treat <= 2018 & year < first_treat ///
    if !missing(first_treat)
replace pre_treated = 0 if missing(pre_treated)
gen byte pre_control = first_treat == 0 | first_treat > 2018 if !missing(first_treat)
replace pre_control = 0 if missing(pre_control)

bysort serie co_comp year: egen byte year_has_pre_treated = max(pre_treated)
gen byte sample_pre = year_has_pre_treated == 1 & !missing(medprof_std) & ///
    (pre_treated == 1 | pre_control == 1)

local graphs ""

foreach serie in "9EF" "EM3" {
foreach cocomp in 1 2 {
    if `cocomp' == 1 local disc "LP"
    if `cocomp' == 2 local disc "Mat"
    if `cocomp' == 1 local titdisc "Portuguese Language"
    if `cocomp' == 2 local titdisc "Mathematics"
    if "`serie'" == "9EF" local titserie "9th Grade"
    if "`serie'" == "EM3" local titserie "3rd Grade of High School"

    qui count if serie == "`serie'" & co_comp == `cocomp' & sample_pre == 1 & pre_treated == 1
    local n_trat = r(N)
    qui count if serie == "`serie'" & co_comp == `cocomp' & sample_pre == 1 & pre_control == 1
    local n_ctrl = r(N)

    if `n_trat' > 10 & `n_ctrl' > 10 {
        twoway ///
            (kdensity medprof_std if serie == "`serie'" & co_comp == `cocomp' & ///
                sample_pre == 1 & pre_control == 1, lcolor(navy) lwidth(medthick)) ///
            (kdensity medprof_std if serie == "`serie'" & co_comp == `cocomp' & ///
                sample_pre == 1 & pre_treated == 1, lcolor(maroon) lwidth(medthick)), ///
            title("`titserie' - `titdisc'", size(medsmall)) ///
            xtitle("Standardized proficiency") ///
            ytitle("Density") ///
            legend(order(1 "Control" 2 "Treated, pre-adoption") ///
                rows(1) ring(0) position(6) size(small) region(style(none))) ///
            graphregion(color(white)) bgcolor(white)
        graph rename dens_pre_`serie'_`disc', replace
        local graphs "`graphs' dens_pre_`serie'_`disc'"
        graph export "$fig/densidade_medprof_pre_`serie'_`disc'.pdf", replace
        copy "$fig/densidade_medprof_pre_`serie'_`disc'.pdf" "$wpfig/densidade_medprof_pre_`serie'_`disc'.pdf", replace
    }
}
}

if "`graphs'" != "" {
    graph combine `graphs', cols(2) ///
        title("Density of proficiency in the pre-treatment period") ///
        graphregion(color(white))
    graph export "$fig/densidade_medprof_pre_trat_controle.pdf", replace
    copy "$fig/densidade_medprof_pre_trat_controle.pdf" "$wpfig/densidade_medprof_pre_trat_controle.pdf", replace
}
