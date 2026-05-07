/*===========================================================================
  06b_export_appendix_figures_english.do
  Re-export appendix figures only, using existing appendix model results.
===========================================================================*/

clear all
set more off

global root     "d:/Andreas/Teste Claude/TCC"
global out      "$root/checkpoints"
global fig      "$root/figuras"
global appfig   "$root/TCC___Andreas_Azambuja_Barbisan/tabelas/imagens/apendice"
global wpappfig "$root/WP_Saresp_Sant_Anna___Xu/files/images/appendix"

capture mkdir "$appfig"
capture mkdir "$wpappfig"

copy "$fig/eventstudy_all_9EF_LP.pdf"  "$appfig/eventstudy_base_9EF_LP.pdf", replace
copy "$fig/eventstudy_all_9EF_LP.pdf"  "$wpappfig/eventstudy_base_9EF_LP.pdf", replace
copy "$fig/eventstudy_all_9EF_Mat.pdf" "$appfig/eventstudy_base_9EF_Mat.pdf", replace
copy "$fig/eventstudy_all_9EF_Mat.pdf" "$wpappfig/eventstudy_base_9EF_Mat.pdf", replace
copy "$fig/eventstudy_all_EM3_LP.pdf"  "$appfig/eventstudy_base_EM3_LP.pdf", replace
copy "$fig/eventstudy_all_EM3_LP.pdf"  "$wpappfig/eventstudy_base_EM3_LP.pdf", replace
copy "$fig/eventstudy_all_EM3_Mat.pdf" "$appfig/eventstudy_base_EM3_Mat.pdf", replace
copy "$fig/eventstudy_all_EM3_Mat.pdf" "$wpappfig/eventstudy_base_EM3_Mat.pdf", replace

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
            legend(order(2 "Main specification" 4 "P10-P90 support") ///
                rows(1) ring(0) position(6) size(small) region(style(none))) ///
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
                legend(order(2 "Main specification" 4 "No imputation") ///
                    rows(1) ring(0) position(6) size(small) region(style(none))) ///
                graphregion(color(white)) bgcolor(white)
            graph export "$appfig/eventstudy_imputacao_`serie'_`disc'.pdf", replace
            copy "$appfig/eventstudy_imputacao_`serie'_`disc'.pdf" "$wpappfig/eventstudy_imputacao_`serie'_`disc'.pdf", replace
        }
    }
}
