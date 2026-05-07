/*===========================================================================
  02_estimacao.do - TWFE, CS21, DR-DiD, S&X aproximado (Nivel B)
  Andreas Azambuja Barbisan - TCC Insper 2026

  VERSAO INTEGRADA E CORRIGIDA (abril 2026)

  Sequencia de estimadores:
    1. TWFE   - benchmark sob adocao escalonada
    2. CS21   - Callaway & Sant'Anna (2021), sem covariadas
    3. DRDID  - Sant'Anna & Zhao (2020), covariadas baseline + composicao
    4. DRRAW  - robustez: covariadas brutas sem imputacao
    5. SXapr  - Sant'Anna & Xu (2026), implementacao parametrica Nivel B

  Sobre o Nivel B (bloco 5):
    - PS multinomial 4-cat (mlogit), nao logit binario
    - Tres OR regressions separadas: m_{1,0}, m_{0,1}, m_{0,0}
    - Estimando tau_dr via eq. (2.6) do artigo, nao via 'drdid'
    - Covariadas fixadas em g-1 da coorte (evita bad controls)
    - Variancia agregada via funcao de influencia empirica, cluster em codesc
    - Sem fallback binario: estimativas S&X reportadas sempre usam PS 4-cat
    - Escada adaptativa de covariadas padronizadas para lidar com coortes raras
    - Verbose: imprime status de cada par (g, l) para diagnostico
===========================================================================*/

clear all
set more off

global root   "d:/Andreas/Teste Claude/TCC"
global saresp "$root/Agregado Saresp"
global out    "$root/checkpoints"
global fig    "$root/figuras"
capture mkdir "$fig"
capture mkdir "$out/logs"

* Janela do event study
global K_PRE  4
global K_POST 4

* Covariadas
global cov_mun  "ln_pib_pc ln_pop prop_jovem_0_19"
global cov_base "medprof_base $cov_mun"
global cov_comp "pct_fem_imp pct_nee_imp idade_media_imp pct_preta_parda_imp pct_raca_declarada_imp pct_bolsa_fam_imp ln_alunos"
global cov_sx   "$cov_comp"
global cov_raw  "pct_fem pct_nee idade_media pct_preta_parda pct_raca_declarada pct_bolsa_fam ln_alunos_raw"

capture log close _all
local runstamp = subinstr("`c(current_date)'_`c(current_time)'", " ", "_", .)
local runstamp = subinstr("`runstamp'", ":", "", .)
global estim_log "$out/logs/02_estimacao_console_`runstamp'.log"
log using "$estim_log", text replace
di "Log desta rodada: $estim_log"

local event_vars ""
forvalues k = $K_PRE(-1)1 {
    local event_vars "`event_vars' Tm`k'"
}
forvalues k = 0/$K_POST {
    local event_vars "`event_vars' Tp`k'"
}
local n_event : word count `event_vars'

*===========================================================================
* PROGRAMAS AUXILIARES
*===========================================================================

capture program drop build_event_est
program define build_event_est, eclass
    syntax name(name=src), Store(name) Eventvars(string) [Nobs(string) Nclust(string)]

    quietly estimates restore `src'

    tempname bnew Vnew

    if "`nobs'" == "" {
        capture local nobs = e(N)
    }
    if "`nclust'" == "" {
        capture local nclust = e(N_clust)
    }

    local kept_eventvars ""
    foreach coef of local eventvars {
        capture local btest = _b[`coef']
        if !_rc {
            local kept_eventvars "`kept_eventvars' `coef'"
        }
        else if "`coef'" == "Tm1" {
            local kept_eventvars "`kept_eventvars' `coef'"
        }
    }

    local k : word count `kept_eventvars'
    matrix `bnew' = J(1, `k', .)
    matrix colnames `bnew' = `kept_eventvars'
    matrix `Vnew' = J(`k', `k', 0)
    matrix colnames `Vnew' = `kept_eventvars'
    matrix rownames `Vnew' = `kept_eventvars'

    local j = 1
    foreach coef of local kept_eventvars {
        capture local bval = _b[`coef']
        if _rc {
            if "`coef'" == "Tm1" {
                local bval = 0
                local sval = 0
            }
            else {
                local bval = .
                local sval = .
            }
        }
        else {
            capture local sval = _se[`coef']
            if _rc local sval = .
        }

        matrix `bnew'[1, `j'] = `bval'
        matrix `Vnew'[`j', `j'] = (`sval')^2
        local ++j
    }

    ereturn post `bnew' `Vnew'
    ereturn local cmd "build_event_est"
    capture confirm number `nobs'
    if !_rc ereturn scalar N = `nobs'
    capture confirm number `nclust'
    if !_rc ereturn scalar N_clust = `nclust'
    estimates store `store'
end

capture program drop build_sx_est
program define build_sx_est, eclass
    syntax , Bmat(name) Vmat(name) Store(name) [Nobs(string) Nclust(string) Eventcells(string)]

    * Se todos os elementos de bmat forem missing, nao post - apenas gera matriz vazia.
    * Isso evita r(504) e permite que post_event_rows depois lide com "sem estimativa".
    local k_cols = colsof(`bmat')
    local n_valid = 0
    forvalues j = 1/`k_cols' {
        if !missing(`bmat'[1, `j']) local ++n_valid
    }

    * Pega os nomes das colunas da matriz original para preservar no output
    local cnames : colnames `bmat'

    if `n_valid' == 0 {
        * Posta um vetor de zeros com variancia zero, so para manter o store.
        tempname bz Vz
        matrix `bz' = J(1, `k_cols', 0)
        matrix colnames `bz' = `cnames'
        matrix `Vz' = J(`k_cols', `k_cols', 0)
        matrix colnames `Vz' = `cnames'
        matrix rownames `Vz' = `cnames'
        ereturn post `bz' `Vz'
        ereturn local cmd "sx_event_empty"
    }
    else {
        * Substitui missings por 0 na diagonal de V e 0 em b para evitar r(504).
        tempname bfix Vfix
        matrix `bfix' = `bmat'
        matrix `Vfix' = `vmat'
        forvalues j = 1/`k_cols' {
            if missing(`bfix'[1, `j']) {
                matrix `bfix'[1, `j'] = 0
                matrix `Vfix'[`j', `j'] = 0
            }
        }
        ereturn post `bfix' `Vfix'
        ereturn local cmd "sx_event"
    }

    capture confirm number `nobs'
    if !_rc ereturn scalar N = `nobs'
    capture confirm number `nclust'
    if !_rc ereturn scalar N_clust = `nclust'
    capture confirm number `eventcells'
    if !_rc ereturn scalar N_event_cells = `eventcells'
    estimates store `store'
end

capture program drop post_event_rows
program define post_event_rows
    syntax name(name=src), Handle(name) Serie(string) Disc(string) Estimator(string) Eventvars(string)

    quietly estimates restore `src'
    capture local nobs = e(N)
    if _rc local nobs = .
    capture local nclust = e(N_clust)
    if _rc local nclust = .
    capture local eventcells = e(N_event_cells)
    if _rc local eventcells = .

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

        local ci_low  = .
        local ci_high = .
        if !missing(`att') & !missing(`se') & `se' > 0 {
            local ci_low  = `att' - 1.96 * `se'
            local ci_high = `att' + 1.96 * `se'
        }

        post `handle' ("`serie'") ("`disc'") ("`estimator'") ///
            (`event_time') (`att') (`se') (`ci_low') (`ci_high') (`nobs') (`nclust') (`eventcells')
    }
end

*===========================================================================
* PREPARAR AMOSTRA DE ESTIMACAO
*===========================================================================

use "$out/painel_completo.dta", clear

* Padronizar proficiencia por ano: z-score dentro de cada ano
capture drop medprof_std medprof_year_mean medprof_year_sd
bysort year: egen double medprof_year_mean = mean(medprof)
bysort year: egen double medprof_year_sd   = sd(medprof)
gen double medprof_std = (medprof - medprof_year_mean) / medprof_year_sd ///
    if !missing(medprof) & medprof_year_sd > 0

* Escolas tratadas apos 2018 entram como not-yet-treated na janela 2011-2018
replace gvar    = 0 if first_treat > 2018 & !missing(first_treat)
replace treated = 0 if first_treat > 2018 & !missing(first_treat)

tempname post_event post_diag post_sx
postfile `post_event' str6 serie str4 disc str8 estimator int event_time ///
    double att se ci_low ci_high N_obs N_clust N_event_cells ///
    using "$out/estimacao_eventstudy_long.dta", replace
postfile `post_diag' str6 serie str4 disc str8 estimator ///
    double N_obs n_schools treated_obs control_obs treated_schools n_cohorts N_clust N_event_cells ///
    str120 notes using "$out/estimacao_diagnosticos.dta", replace
postfile `post_sx' str6 serie str4 disc int rel_time int cohort int cal_time ///
    double n_sample n_treated n11 n10 n01 n00 byte converged double att se ///
    double p11_min p11_p1 p10_p1 p01_p1 p00_p1 cut11 cut10 cut01 cut00 ///
    str32 status using "$out/estimacao_sx_detalhe.dta", replace

*===========================================================================
* LOOP PRINCIPAL: serie x disciplina
*===========================================================================

foreach serie in "9EF" "EM3" {
foreach cocomp in 1 2 {

    if `cocomp' == 1 local disc "LP"
    if `cocomp' == 2 local disc "Mat"

    di _newline(2) "{hline 70}"
    di "SERIE: `serie' | DISCIPLINA: `disc'"
    di "{hline 70}"

    preserve
    keep if serie == "`serie'" & co_comp == `cocomp'

    tempvar school_tag treated_school
    bysort codesc: gen byte `school_tag' = (_n == 1)
    bysort codesc: egen byte `treated_school' = max(treated == 1)

    qui count
    local N_cell = r(N)
    qui count if `school_tag'
    local n_schools = r(N)
    qui count if treated == 1
    local treated_obs = r(N)
    qui count if treated == 0
    local control_obs = r(N)
    qui count if `school_tag' & `treated_school' == 1
    local treated_schools = r(N)

    di "  Amostra: N=" %10.0f `N_cell' " | escolas=" %7.0f `n_schools' ///
       " | obs tratadas=" %10.0f `treated_obs' " | obs controle=" %10.0f `control_obs'

    qui count if treated == 1 & !missing(medprof_std)
    if r(N) < 30 {
        di "  Poucos tratados com medprof_std observada (" %9.0f r(N) "), pulando."
        post `post_diag' ("`serie'") ("`disc'") ("SKIP") (`N_cell') (`n_schools') ///
            (`treated_obs') (`control_obs') (`treated_schools') (.) (.) (.) ///
            ("menos de 30 tratados com medprof_std")
        restore
        continue
    }

    xtset codesc year

    *-------------------------------------------------------------------
    * 1. TWFE
    *-------------------------------------------------------------------
    di _newline "--- 1. TWFE ---"

    gen rel_twfe = max(-$K_PRE, min($K_POST, rel_time)) if treated == 1
    replace rel_twfe = -999 if treated == 0 | missing(rel_time)

    forvalues k = $K_PRE(-1)2 {
        gen Tm`k' = (rel_twfe == -`k')
    }
    forvalues k = 0/$K_POST {
        gen Tp`k' = (rel_twfe == `k')
    }

    local pre_dummies ""
    forvalues k = $K_PRE(-1)2 {
        local pre_dummies "`pre_dummies' Tm`k'"
    }
    local post_dummies ""
    forvalues k = 0/$K_POST {
        local post_dummies "`post_dummies' Tp`k'"
    }

    reghdfe medprof_std `pre_dummies' `post_dummies' $cov_sx, ///
        absorb(codesc year) cluster(codesc)
    estimates store twfe_`serie'_`disc'
    local twfe_N     = e(N)
    local twfe_clust = e(N_clust)

    build_event_est twfe_`serie'_`disc', ///
        store(twfe_plot_`serie'_`disc') eventvars("`event_vars'") ///
        nobs(`twfe_N') nclust(`twfe_clust')

    post `post_diag' ("`serie'") ("`disc'") ("TWFE") (`twfe_N') (`n_schools') ///
        (`treated_obs') (`control_obs') (`treated_schools') (.) (`twfe_clust') (.) ///
        ("FE escola+ano; cov_comp em t")

    drop rel_twfe Tm* Tp*

    *-------------------------------------------------------------------
    * 2. CS21 - sem covariadas
    *-------------------------------------------------------------------
    di _newline "--- 2. CS21 (sem covariadas) ---"

    csdid medprof_std, ivar(codesc) time(year) gvar(gvar) ///
        method(dripw) notyet long2
    local cs_N = e(N)
    local cs_clust = .
    capture local cs_clust = e(N_clust)
    estat event, window(-$K_PRE $K_POST) estore(cs_nc_`serie'_`disc')

    estimates restore cs_nc_`serie'_`disc'

    build_event_est cs_nc_`serie'_`disc', ///
        store(cs_nc_plot_`serie'_`disc') eventvars("`event_vars'") ///
        nobs(`cs_N') nclust(`cs_clust')

    *-------------------------------------------------------------------
    * 3. DR-DiD - covariadas baseline + composicao agregada
    *-------------------------------------------------------------------
    di _newline "--- 3. DR-DiD (covariadas baseline) ---"

    csdid medprof_std $cov_base $cov_sx, ivar(codesc) time(year) gvar(gvar) ///
        method(dripw) notyet long2
    local drdid_N = e(N)
    local drdid_clust = .
    capture local drdid_clust = e(N_clust)
    estat event, window(-$K_PRE $K_POST) estore(drdid_`serie'_`disc')

    estimates restore drdid_`serie'_`disc'

    build_event_est drdid_`serie'_`disc', ///
        store(drdid_plot_`serie'_`disc') eventvars("`event_vars'") ///
        nobs(`drdid_N') nclust(`drdid_clust')

    *-------------------------------------------------------------------
    * 3b. Robustez - DR-DiD com covariadas brutas, sem imputacao
    *-------------------------------------------------------------------
    di _newline "--- 3b. DR-DiD robustez (sem imputacao) ---"

    local drraw_N = .
    local drraw_clust = .
    capture noisily csdid medprof_std $cov_base $cov_raw, ivar(codesc) time(year) gvar(gvar) ///
        method(dripw) notyet long2
    local rc_drraw = _rc
    if !`rc_drraw' {
        local drraw_N = e(N)
        capture local drraw_clust = e(N_clust)
        capture noisily estat event, window(-$K_PRE $K_POST) estore(drraw_`serie'_`disc')
        local rc_drraw = _rc
        if !`rc_drraw' {
            estimates restore drraw_`serie'_`disc'
            build_event_est drraw_`serie'_`disc', ///
                store(drraw_plot_`serie'_`disc') eventvars("`event_vars'") ///
                nobs(`drraw_N') nclust(`drraw_clust')
        }
        else {
            di "  Robustez sem imputacao estimou csdid, mas estat event falhou: rc=`rc_drraw'"
        }
    }
    else {
        di "  Robustez sem imputacao falhou: rc=`rc_drraw'"
    }

    qui levelsof gvar if treated == 1 & gvar <= 2017, local(cohorts)
    local n_cohorts : word count `cohorts'

    post `post_diag' ("`serie'") ("`disc'") ("CS21") (`cs_N') (`n_schools') ///
        (`treated_obs') (`control_obs') (`treated_schools') (`n_cohorts') (`cs_clust') (.) ///
        ("csdid dripw notyet sem covariadas")
    post `post_diag' ("`serie'") ("`disc'") ("DRDID") (`drdid_N') (`n_schools') ///
        (`treated_obs') (`control_obs') (`treated_schools') (`n_cohorts') (`drdid_clust') (.) ///
        ("csdid dripw notyet com cov_base e cov_comp")
    post `post_diag' ("`serie'") ("`disc'") ("DRRAW") (`drraw_N') (`n_schools') ///
        (`treated_obs') (`control_obs') (`treated_schools') (`n_cohorts') (`drraw_clust') (.) ///
        ("robustez: covariadas brutas sem imputacao; rc=`rc_drraw'")

    *-------------------------------------------------------------------
    * 4. S&X aproximado (Nivel B): mlogit 4-cat + 3 OR separadas
    *-------------------------------------------------------------------
    di _newline "--- 4. S&X aproximado (Nivel B) ---"
    di "  Coortes ativas: `cohorts'"

    local lmin = -$K_PRE
    local lmax =  $K_POST
    local nL   = $K_PRE + $K_POST + 1

    matrix sx_l   = J(`nL', 1, .)
    matrix sx_att = J(`nL', 1, .)
    matrix sx_se  = J(`nL', 1, .)

    local sx_pairs_ok_total = 0
    local sx_rc_fail_total  = 0
    local sx_lowobs_total   = 0
    local sx_mlogit_fail    = 0
    local sx_ps_fail        = 0
    local sx_reduced_used   = 0

    * Covariadas para o PS multinomial e as OR.
    * A estimacao S&X abaixo nao usa fallback binario. Se o mlogit 4-cat
    * rico nao converge, tentamos uma escada parcimoniosa ainda multinomial.
    local cov_sx_list "pct_fem_imp idade_media_imp pct_preta_parda_imp pct_bolsa_fam_imp ln_alunos"
    local k_cov : word count `cov_sx_list'
    local sx_model_1 "pct_fem_imp idade_media_imp pct_preta_parda_imp pct_bolsa_fam_imp ln_alunos"
    local sx_model_2 "idade_media_imp pct_preta_parda_imp pct_bolsa_fam_imp ln_alunos"
    local sx_model_3 "idade_media_imp pct_bolsa_fam_imp ln_alunos"
    local sx_model_4 "pct_bolsa_fam_imp ln_alunos"
    local sx_model_5 "ln_alunos"

    local row = 1
    forvalues l = `lmin'(1)`lmax' {
        matrix sx_l[`row', 1] = `l'

        if `l' == -1 {
            di "      l=" %3.0f `l' " -> baseline universal normalizado para zero"
            foreach g of local cohorts {
                local t = `g' + `l'
                if `t' < 2011 | `t' > 2018 continue
                post `post_sx' ("`serie'") ("`disc'") (`l') (`g') (`t') ///
                    (.) (.) (.) (.) (.) (.) (0) (.) (.) ///
                    (.) (.) (.) (.) (.) (.) (.) (.) (.) ///
                    ("baseline_km1")
            }
            local row = `row' + 1
            continue
        }

        local num_w_att = 0
        local den_w     = 0
        local var_agg   = 0
        local sx_pairs_ok = 0

        foreach g of local cohorts {
            local t = `g' + `l'
            if `t' < 2011 | `t' > 2018 continue

            qui count if (gvar == `g' | gvar == 0) ///
                & inlist(year, `g' - 1, `t') & !missing(medprof_std)
            local n_pair_raw = r(N)

            * Precisamos de n suficiente para mlogit de 4-cat com k_cov covariadas.
            * Heuristica: pelo menos 15 + 3*k_cov observacoes por par.
            local n_min_total = 15 + 3 * `k_cov'
            if `n_pair_raw' < `n_min_total' {
                local ++sx_lowobs_total
                post `post_sx' ("`serie'") ("`disc'") (`l') (`g') (`t') ///
                    (`n_pair_raw') (.) (.) (.) (.) (.) (0) (.) (.) ///
                    (.) (.) (.) (.) (.) (.) (.) (.) (.) ///
                    ("n_too_small")
                continue
            }

            local att_g = .
            local se_g  = .
            local n_11  = 0
            local n_10  = 0
            local n_01  = 0
            local n_00  = 0
            local n_pair = 0
            local converged = 0
            local status "init"
            local p11_min = .
            local p11_p1  = .
            local p10_p1  = .
            local p01_p1  = .
                local p00_p1  = .
                local cut11   = .
                local cut10   = .
                local cut01   = .
                local cut00   = .
                local sx_model_used ""
                local sx_k_used = .

            capture frame drop gt_frame
            frame copy default gt_frame
            frame change gt_frame

                quietly keep if (gvar == `g' | gvar == 0) & inlist(year, `g' - 1, `t')

                * Variaveis do setup 2x2 (criadas antes do baseline para uso nos egens)
                quietly gen byte D_sx    = (gvar == `g')
                quietly gen byte period  = (year == `t')
                quietly gen byte cat_11  = (D_sx == 1) & (period == 1)
                quietly gen byte cat_10  = (D_sx == 1) & (period == 0)
                quietly gen byte cat_01  = (D_sx == 0) & (period == 1)
                quietly gen byte cat_00  = (D_sx == 0) & (period == 0)

                * Variavel categorica para mlogit: 0=(0,0), 1=(0,1), 2=(1,0), 3=(1,1)
                quietly gen byte dtcat = .
                quietly replace dtcat = 0 if cat_00 == 1
                quietly replace dtcat = 1 if cat_01 == 1
                quietly replace dtcat = 2 if cat_10 == 1
                quietly replace dtcat = 3 if cat_11 == 1

                * --- BASELINE COVARIATES: fixa covariadas em g-1 ---
                * Para cada escola, pega valores das covariadas no ano g-1
                * e replica para o ano t. Evita bad controls.
                local cov_base_list ""
                foreach v of local cov_sx_list {
                    quietly gen double `v'_b = `v' if year == `g' - 1
                    quietly bysort codesc (year): egen double `v'_base = max(`v'_b)
                    quietly drop `v'_b
                    local cov_base_list "`cov_base_list' `v'_base"
                }

                * Dropar observacoes sem baseline (escolas que so aparecem no pos)
                * ou sem medprof_std observada
                quietly egen byte miss_any = rowmiss(medprof_std `cov_base_list' dtcat)
                quietly keep if miss_any == 0
                quietly drop miss_any

                qui count
                local n_pair = r(N)
                qui count if dtcat == 3
                local n_11 = r(N)
                qui count if dtcat == 2
                local n_10 = r(N)
                qui count if dtcat == 1
                local n_01 = r(N)
                qui count if dtcat == 0
                local n_00 = r(N)

                * Minimo por celula. A escada adaptativa pode usar modelos
                * menores quando as coortes do EM3 sao pequenas.
                local n_min_cell = 5
                if `n_11' < `n_min_cell' | `n_10' < `n_min_cell' | ///
                   `n_01' < `n_min_cell' | `n_00' < `n_min_cell' {
                    local status "cell_too_small"
                    local att_g = .
                    local se_g  = .
                }
                else {
                    * Padronizar covariadas dentro do par 2x2. Isso ajuda o
                    * otimizador do mlogit e evita que escala numerica vire
                    * falsa nao-convergencia.
                    foreach v of local cov_sx_list {
                        quietly summarize `v'_base
                        local mu = r(mean)
                        local sd = r(sd)
                        if missing(`sd') | `sd' <= 1e-8 {
                            quietly gen double z_`v'_base = 0
                        }
                        else {
                            quietly gen double z_`v'_base = (`v'_base - `mu') / `sd'
                        }
                    }

                    * ---------------------------------------------------------
                    * (a) PS MULTINOMIAL: mlogit 4-cat, base=(0,0).
                    *     Sem fallback para S&Z: se entrar, veio de PS 4-cat.
                    * ---------------------------------------------------------
                    local conv_ml = 0
                    local cov_model ""

                    forvalues sx_try = 1/5 {
                        if `conv_ml' == 0 {
                            local cov_try ""
                            foreach v of local sx_model_`sx_try' {
                                local cov_try "`cov_try' z_`v'_base"
                            }

                            capture quietly mlogit dtcat `cov_try', base(0) ///
                                difficult technique(nr bhhh dfp) iterate(300)
                            local rc_ml = _rc
                            local conv_try = 0
                            if !`rc_ml' capture local conv_try = e(converged)

                            if !`rc_ml' & `conv_try' == 1 {
                                local conv_ml = 1
                                local cov_model "`cov_try'"
                                local sx_k_used : word count `cov_try'
                                local sx_model_used "k`sx_k_used'"
                                if `sx_try' > 1 local ++sx_reduced_used
                            }
                        }
                    }

                    if `conv_ml' == 1 {
                        quietly predict double p00_hat, outcome(0)
                        quietly predict double p01_hat, outcome(1)
                        quietly predict double p10_hat, outcome(2)
                        quietly predict double p11_hat, outcome(3)

                        * Checagem de overlap focada: 
                        * A regra original com minimo absoluto estava dura demais
                        * para o 9EF, onde a celula (1,1) e rara. Passamos a usar
                        * o percentil 1 e cortes relativos ao tamanho observado de
                        * cada celula, agora sem fallback binario.
                        qui summarize p11_hat, meanonly
                        local p11_min = r(min)
                        capture quietly _pctile p11_hat, p(1)
                        if !_rc local p11_p1 = r(r1)
                        if missing(`p11_p1') local p11_p1 = `p11_min'

                        local share11 = `n_11' / `n_pair'
                        local share10 = `n_10' / `n_pair'
                        local share01 = `n_01' / `n_pair'
                        local share00 = `n_00' / `n_pair'

                        local cut11 = max(0.0005, 0.10 * `share11')
                        local cut10 = max(0.0005, 0.10 * `share10')
                        local cut01 = max(0.0005, 0.10 * `share01')
                        local cut00 = max(0.0005, 0.10 * `share00')

                        local ps_ok = (`p11_p1' >= `cut11')
                        if `ps_ok' {
                            qui summarize p10_hat if cat_10 == 1, meanonly
                            local n_ps10 = r(N)
                            local p10_min = r(min)
                            capture quietly _pctile p10_hat if cat_10 == 1, p(1)
                            if !_rc local p10_p1 = r(r1)
                            if missing(`p10_p1') local p10_p1 = `p10_min'
                            if `n_ps10' > 0 local ps_ok = `ps_ok' & (`p10_p1' >= `cut10')
                        }
                        if `ps_ok' {
                            qui summarize p01_hat if cat_01 == 1, meanonly
                            local n_ps01 = r(N)
                            local p01_min = r(min)
                            capture quietly _pctile p01_hat if cat_01 == 1, p(1)
                            if !_rc local p01_p1 = r(r1)
                            if missing(`p01_p1') local p01_p1 = `p01_min'
                            if `n_ps01' > 0 local ps_ok = `ps_ok' & (`p01_p1' >= `cut01')
                        }
                        if `ps_ok' {
                            qui summarize p00_hat if cat_00 == 1, meanonly
                            local n_ps00 = r(N)
                            local p00_min = r(min)
                            capture quietly _pctile p00_hat if cat_00 == 1, p(1)
                            if !_rc local p00_p1 = r(r1)
                            if missing(`p00_p1') local p00_p1 = `p00_min'
                            if `n_ps00' > 0 local ps_ok = `ps_ok' & (`p00_p1' >= `cut00')
                        }

                        if `ps_ok' {
                            * ---------------------------------------------------------
                            * (b) OR REGRESSIONS: tres regressoes separadas
                            * ---------------------------------------------------------
                            capture quietly regress medprof_std `cov_model' if cat_10 == 1
                            local rc_or10 = _rc
                            if !`rc_or10' quietly predict double m10_hat, xb

                            capture quietly regress medprof_std `cov_model' if cat_01 == 1
                            local rc_or01 = _rc
                            if !`rc_or01' quietly predict double m01_hat, xb

                            capture quietly regress medprof_std `cov_model' if cat_00 == 1
                            local rc_or00 = _rc
                            if !`rc_or00' quietly predict double m00_hat, xb

                            if !`rc_or10' & !`rc_or01' & !`rc_or00' {
                                * ---------------------------------------------------------
                                * (c) PESOS HAJEK (eq. 2.4 do artigo)
                                * ---------------------------------------------------------
                                quietly gen double w11_raw = cat_11
                                quietly gen double w10_raw = cat_10 * p11_hat / p10_hat
                                quietly gen double w01_raw = cat_01 * p11_hat / p01_hat
                                quietly gen double w00_raw = cat_00 * p11_hat / p00_hat

                                qui summarize w11_raw, meanonly
                                local E_w11 = r(mean)
                                qui summarize w10_raw, meanonly
                                local E_w10 = r(mean)
                                qui summarize w01_raw, meanonly
                                local E_w01 = r(mean)
                                qui summarize w00_raw, meanonly
                                local E_w00 = r(mean)

                                if `E_w11' > 0 & `E_w10' > 0 & `E_w01' > 0 & `E_w00' > 0 {
                                    quietly gen double w11 = w11_raw / `E_w11'
                                    quietly gen double w10 = w10_raw / `E_w10'
                                    quietly gen double w01 = w01_raw / `E_w01'
                                    quietly gen double w00 = w00_raw / `E_w00'

                                    * ---------------------------------------------------------
                                    * (d) ESTIMAND tau_dr via eq. (2.6) do artigo
                                    *
                                    * Sinais: (-1)^(d+t) dentro do somatorio -
                                    *   (1,0): -1, (0,1): -1, (0,0): +1
                                    * ---------------------------------------------------------
                                    quietly gen double psi_i = ///
                                        w11 * (medprof_std - m10_hat - m01_hat + m00_hat) ///
                                      - w10 * (medprof_std - m10_hat) ///
                                      - w01 * (medprof_std - m01_hat) ///
                                      + w00 * (medprof_std - m00_hat)

                                    qui summarize psi_i, meanonly
                                    local att_g = r(mean)

                                    * ---------------------------------------------------------
                                    * (e) VARIANCIA clusterizada em codesc
                                    *     V_g = (1/n^2) * sum_{g in clusters} (sum_{i in g} psi_c_i)^2
                                    * ---------------------------------------------------------
                                    quietly gen double psi_c = psi_i - `att_g'
                                    quietly bysort codesc: egen double psi_cluster = total(psi_c)
                                    quietly bysort codesc: gen byte school_tag2 = (_n == 1)
                                    qui count
                                    local n_total = r(N)
                                    qui summarize psi_cluster if school_tag2
                                    local sum_sq = 0
                                    qui gen double psi_cluster_sq = psi_cluster^2 if school_tag2
                                    qui summarize psi_cluster_sq, meanonly
                                    local sum_sq = r(sum)
                                    if `n_total' > 0 {
                                        local var_g = `sum_sq' / (`n_total'^2)
                                        local se_g  = sqrt(`var_g')
                                        local converged = 1
                                        local status "ok_mlogit_`sx_model_used'"
                                    }
                                    else {
                                        local status "var_fail"
                                        local att_g = .
                                        local se_g  = .
                                    }
                                }
                                else {
                                    local status "hajek_zero"
                                    local att_g = .
                                    local se_g  = .
                                }
                            }
                            else {
                                local status "or_fail"
                                local ++sx_rc_fail_total
                                local att_g = .
                                local se_g  = .
                            }
                        }
                        else {
                            local ++sx_ps_fail
                            local status "ps_lowoverlap_`sx_model_used'"
                            local att_g = .
                            local se_g  = .
                        }
                    }
                    else {
                        * mlogit nao convergiu ou falhou. Mantemos ausente em
                        * vez de substituir por S&Z binario, para preservar a
                        * interpretacao S&X do estimador reportado.
                        local ++sx_mlogit_fail
                        local status "mlogit_fail_no_sz"
                        local att_g = .
                        local se_g  = .
                    }
                }

            frame change default
            frame drop gt_frame

            * Impressao verbose por coorte
            if !missing(`att_g') {
                di "      g=" %4.0f `g' " l=" %3.0f `l' " -> att=" %7.3f `att_g' ///
                   " se=" %6.3f `se_g' " | n11=" %4.0f `n_11' " n10=" %4.0f `n_10' ///
                   " n01=" %4.0f `n_01' " n00=" %4.0f `n_00' ///
                   " | p11_q1=" %6.4f `p11_p1' " cut11=" %6.4f `cut11' ///
                   " [" "`status'" "]"
            }
            else {
                di "      g=" %4.0f `g' " l=" %3.0f `l' " -> FAIL" ///
                   " | n11=" %4.0f `n_11' " n10=" %4.0f `n_10' ///
                   " n01=" %4.0f `n_01' " n00=" %4.0f `n_00' ///
                   " | p11_q1=" %6.4f `p11_p1' " cut11=" %6.4f `cut11' ///
                   " [" "`status'" "]"
            }

            * Agregacao tipo C&S: peso = n_11 (tratados no pos da coorte)
            if !missing(`att_g') & !missing(`se_g') & `n_11' > 0 {
                local num_w_att = `num_w_att' + `att_g' * `n_11'
                local den_w     = `den_w'     + `n_11'
                local var_agg   = `var_agg' + (`n_11')^2 * (`se_g')^2
                local ++sx_pairs_ok
                local ++sx_pairs_ok_total
            }

            post `post_sx' ("`serie'") ("`disc'") (`l') (`g') (`t') ///
                (`n_pair') (`n_11') (`n_11') (`n_10') (`n_01') (`n_00') ///
                (`converged') (`att_g') (`se_g') ///
                (`p11_min') (`p11_p1') (`p10_p1') (`p01_p1') (`p00_p1') ///
                (`cut11') (`cut10') (`cut01') (`cut00') ///
                ("`status'")
        }

        if `den_w' > 0 {
            matrix sx_att[`row', 1] = `num_w_att' / `den_w'
            matrix sx_se[`row', 1]  = sqrt(`var_agg') / `den_w'
            di "    => l=" %3.0f `l' " agregado: ATT=" %9.3f (`num_w_att' / `den_w') ///
               " SE=" %9.3f (sqrt(`var_agg') / `den_w') ///
               " coortes validas=" %3.0f `sx_pairs_ok'
        }
        else {
            di "    => l=" %3.0f `l' " sem estimativa agregada valida"
        }
        local row = `row' + 1
    }

    di _newline "  [SXapr] resumo: pairs_ok=`sx_pairs_ok_total', ps_fail=`sx_ps_fail', mlogit_fail=`sx_mlogit_fail', reduced_mlogit=`sx_reduced_used', low_n=`sx_lowobs_total', or_fail=`sx_rc_fail_total'"

    *-------------------------------------------------------------------
    * Salvar resultados S&X aproximado
    *-------------------------------------------------------------------

    capture frame drop sx_frame
    frame create sx_frame
    frame sx_frame {
        clear
        quietly svmat sx_l,   names(l)
        quietly svmat sx_att, names(att)
        quietly svmat sx_se,  names(se)
        quietly rename l1   l
        quietly rename att1 att
        quietly rename se1  se
        quietly gen ci_low  = att - 1.96 * se
        quietly gen ci_high = att + 1.96 * se
        quietly gen serie_disc = "`serie'_`disc'"
        quietly drop if missing(att)
        save "$out/sx_`serie'_`disc'.dta", replace
    }
    frame drop sx_frame

    * Construir matriz b e V no formato do event study
    tempname b_sx V_sx
    matrix `b_sx' = J(1, `n_event', .)
    matrix colnames `b_sx' = `event_vars'
    matrix `V_sx' = J(`n_event', `n_event', 0)
    matrix colnames `V_sx' = `event_vars'
    matrix rownames `V_sx' = `event_vars'

    local j = 1
    foreach coef of local event_vars {
        if "`coef'" == "Tm1" {
            matrix `b_sx'[1, `j'] = 0
            matrix `V_sx'[`j', `j'] = 0
        }
        else {
            if substr("`coef'", 1, 2) == "Tm" {
                local this_l = -real(substr("`coef'", 3, .))
            }
            else {
                local this_l = real(substr("`coef'", 3, .))
            }
            local sx_row = `this_l' - `lmin' + 1
            if `sx_row' >= 1 & `sx_row' <= `nL' {
                local att_now = el(sx_att, `sx_row', 1)
                local se_now  = el(sx_se,  `sx_row', 1)
                if !missing(`att_now') {
                    matrix `b_sx'[1, `j'] = `att_now'
                    if !missing(`se_now') {
                        matrix `V_sx'[`j', `j'] = (`se_now')^2
                    }
                }
            }
        }
        local ++j
    }

    build_sx_est, bmat(`b_sx') vmat(`V_sx') ///
        store(sx_plot_`serie'_`disc') nclust(`n_schools') eventcells(`sx_pairs_ok_total')

    post `post_diag' ("`serie'") ("`disc'") ("SXapr") (.) (`n_schools') ///
        (`treated_obs') (`control_obs') (`treated_schools') (`n_cohorts') ///
        (`n_schools') (`sx_pairs_ok_total') ///
        ("SXapr: mlogit 4-cat + 3 OR; cov baseline; no S&Z fallback; reduced=`sx_reduced_used'")

    *-------------------------------------------------------------------
    * 5. Grafico unico - 4 estimadores no mesmo coefplot
    *-------------------------------------------------------------------
    di _newline "--- 5. Plot e outputs ---"

    if `cocomp' == 1 local titdisc "Língua Portuguesa"
    if `cocomp' == 2 local titdisc "Matemática"
    if "`serie'" == "9EF" local titserie "9º Ano EF"
    if "`serie'" == "EM3" local titserie "3ª Série EM"

    post_event_rows twfe_plot_`serie'_`disc',  handle(`post_event') serie("`serie'") disc("`disc'") estimator("TWFE")   eventvars("`event_vars'")
    post_event_rows cs_nc_plot_`serie'_`disc', handle(`post_event') serie("`serie'") disc("`disc'") estimator("CS21")   eventvars("`event_vars'")
    post_event_rows drdid_plot_`serie'_`disc', handle(`post_event') serie("`serie'") disc("`disc'") estimator("DRDID")  eventvars("`event_vars'")
    capture estimates restore drraw_plot_`serie'_`disc'
    if !_rc {
        post_event_rows drraw_plot_`serie'_`disc', handle(`post_event') serie("`serie'") disc("`disc'") estimator("DRRAW") eventvars("`event_vars'")
    }
    post_event_rows sx_plot_`serie'_`disc',    handle(`post_event') serie("`serie'") disc("`disc'") estimator("SXapr")  eventvars("`event_vars'")

    cap drop zero Tm1
    gen zero = 0
    gen Tm1  = _n
    quietly regress zero Tm1
    estimates store pre_menos_1_`serie'_`disc'
    drop zero Tm1

    di _newline "  Resumo padronizado do event study (baseline comum em k=-1):"
    estimates table twfe_plot_`serie'_`disc' cs_nc_plot_`serie'_`disc' ///
        drdid_plot_`serie'_`disc' sx_plot_`serie'_`disc', ///
        b(%9.3f) se(%9.3f) keep(`event_vars')

    * So faz o coefplot se temos pelo menos uma estimativa S&X valida,
    * caso contrario o drop do Tm1 do sx_plot falha.
    if `sx_pairs_ok_total' > 0 {
        capture noisily coefplot ///
            (twfe_plot_`serie'_`disc'  pre_menos_1_`serie'_`disc', label("TWFE") ///
                msymbol(O) msize(small) offset(-0.18) ///
                lcolor(navy) mcolor(navy) ciopts(recast(rcap) color(navy%90))) ///
            (cs_nc_plot_`serie'_`disc' pre_menos_1_`serie'_`disc', label("CS21") ///
                msymbol(T) msize(small) offset(-0.06) ///
                lcolor(forest_green) mcolor(forest_green) ciopts(recast(rcap) color(forest_green%90))) ///
            (drdid_plot_`serie'_`disc' pre_menos_1_`serie'_`disc', label("DR-DiD") ///
                msymbol(D) msize(small) offset(0.06) ///
                lcolor(maroon) mcolor(maroon) ciopts(recast(rcap) color(maroon%90))) ///
            (sx_plot_`serie'_`disc'    pre_menos_1_`serie'_`disc', label("S&X apr") ///
                msymbol(S) msize(small) offset(0.18) ///
                lcolor(dkorange) mcolor(dkorange) ciopts(recast(rcap) color(dkorange%90))), ///
            vertical yline(0, lcolor(black) lpattern(dash)) ///
            nolabel omitted baselevels keep(*:) ///
            drop(_cons Pre_avg Post_avg ///
                 twfe_plot_`serie'_`disc':Tm1 ///
                 cs_nc_plot_`serie'_`disc':Tm1 ///
                 drdid_plot_`serie'_`disc':Tm1 ///
                 sx_plot_`serie'_`disc':Tm1) ///
            rename(Tm1 = -1 Tm2 = -2 Tm3 = -3 Tm4 = -4 Tp0 = 0 Tp1 = 1 Tp2 = 2 Tp3 = 3 Tp4 = 4) ///
            relocate(-4 = 1 -3 = 2 -2 = 3 -1 = 4 0 = 5 1 = 6 2 = 7 3 = 8 4 = 9) ///
            xline(4.5, lc(gs8) lp(dash)) ///
            title("`titserie' - `titdisc': 4 estimadores") ///
            subtitle("IC 95%") ///
            xtitle("Tempo relativo à conversao PEI (anos)") ///
            ytitle("Efeito médio do tratamento (desvios-padrão)") ///
            legend(order(1 "TWFE" 3 "CS21" 5 "DR-DiD" 7 "S&X apr") ///
                rows(1) ring(0) position(6) xoffset(-7) size(small) region(style(none))) ///
            graphregion(color(white)) bgcolor(white) plotregion(margin(b=8))
        capture graph export "$fig/eventstudy_all_`serie'_`disc'.pdf", replace
        capture graph export "$fig/eventstudy_all_`serie'_`disc'.png", width(2400) replace
    }
    else {
        di "  AVISO: S&X aproximado nao produziu estimativas validas. Gerando coefplot sem S&X."
        capture noisily coefplot ///
            (twfe_plot_`serie'_`disc'  pre_menos_1_`serie'_`disc', label("TWFE") ///
                msymbol(O) msize(small) offset(-0.15) ///
                lcolor(navy) mcolor(navy) ciopts(recast(rcap) color(navy%90))) ///
            (cs_nc_plot_`serie'_`disc' pre_menos_1_`serie'_`disc', label("CS21") ///
                msymbol(T) msize(small) offset(0) ///
                lcolor(forest_green) mcolor(forest_green) ciopts(recast(rcap) color(forest_green%90))) ///
            (drdid_plot_`serie'_`disc' pre_menos_1_`serie'_`disc', label("DR-DiD") ///
                msymbol(D) msize(small) offset(0.15) ///
                lcolor(maroon) mcolor(maroon) ciopts(recast(rcap) color(maroon%90))), ///
            vertical yline(0, lcolor(black) lpattern(dash)) ///
            nolabel omitted baselevels keep(*:) ///
            drop(_cons Pre_avg Post_avg ///
                 twfe_plot_`serie'_`disc':Tm1 ///
                 cs_nc_plot_`serie'_`disc':Tm1 ///
                 drdid_plot_`serie'_`disc':Tm1) ///
            rename(Tm1 = -1 Tm2 = -2 Tm3 = -3 Tm4 = -4 Tp0 = 0 Tp1 = 1 Tp2 = 2 Tp3 = 3 Tp4 = 4) ///
            relocate(-4 = 1 -3 = 2 -2 = 3 -1 = 4 0 = 5 1 = 6 2 = 7 3 = 8 4 = 9) ///
            xline(4.5, lc(gs8) lp(dash)) ///
            title("`titserie' - `titdisc': 3 estimadores (S&X apr falhou)") ///
            subtitle("IC 95%") ///
            xtitle("Tempo relativo à conversao PEI (anos)") ///
            ytitle("Efeito médio do tratamento (desvios-padrão)") ///
            graphregion(color(white)) bgcolor(white)
        capture graph export "$fig/eventstudy_all_`serie'_`disc'.pdf", replace
        capture graph export "$fig/eventstudy_all_`serie'_`disc'.png", width(2400) replace
    }

    restore

} // end co_comp
} // end serie

postclose `post_event'
postclose `post_diag'
postclose `post_sx'

*===========================================================================
* EXPORTS E RESUMOS
*===========================================================================

preserve
    use "$out/estimacao_eventstudy_long.dta", clear
    sort serie disc estimator event_time
    export delimited using "$out/estimacao_eventstudy_long.csv", replace
restore

preserve
    use "$out/estimacao_eventstudy_long.dta", clear

    tempname post_sum
    postfile `post_sum' str6 serie str4 disc str8 estimator ///
        double pre_avg pre_sig_count post_avg k0 k1 k2 k3 k4 N_obs N_event_cells ///
        using "$out/estimacao_resumo_tese.dta", replace

    levelsof serie, local(sum_series)
    foreach s of local sum_series {
        levelsof disc if serie == "`s'", local(sum_discs)
        foreach d of local sum_discs {
            levelsof estimator if serie == "`s'" & disc == "`d'", local(sum_estimators)
            foreach e of local sum_estimators {
                quietly count if serie == "`s'" & disc == "`d'" & estimator == "`e'"
                if r(N) == 0 continue

                quietly summarize att if serie == "`s'" & disc == "`d'" & estimator == "`e'" ///
                    & inlist(event_time, -4, -3, -2), meanonly
                local pre_avg = r(mean)

                quietly count if serie == "`s'" & disc == "`d'" & estimator == "`e'" ///
                    & inlist(event_time, -4, -3, -2) & se > 0 & abs(att / se) > 1.96
                local pre_sig_count = r(N)

                quietly summarize att if serie == "`s'" & disc == "`d'" & estimator == "`e'" ///
                    & event_time >= 0, meanonly
                local post_avg = r(mean)

                foreach kk in 0 1 2 3 4 {
                    quietly summarize att if serie == "`s'" & disc == "`d'" & estimator == "`e'" ///
                        & event_time == `kk', meanonly
                    local k`kk' = r(mean)
                }

                quietly summarize N_obs if serie == "`s'" & disc == "`d'" & estimator == "`e'", meanonly
                local n_obs = r(mean)

                quietly summarize N_event_cells if serie == "`s'" & disc == "`d'" & estimator == "`e'", meanonly
                local n_event_cells = r(mean)

                post `post_sum' ("`s'") ("`d'") ("`e'") ///
                    (`pre_avg') (`pre_sig_count') (`post_avg') ///
                    (`k0') (`k1') (`k2') (`k3') (`k4') (`n_obs') (`n_event_cells')
            }
        }
    }

    postclose `post_sum'
    use "$out/estimacao_resumo_tese.dta", clear
    sort serie disc estimator
    export delimited using "$out/estimacao_resumo_tese.csv", replace
restore

preserve
    use "$out/estimacao_diagnosticos.dta", clear
    sort serie disc estimator
    export delimited using "$out/estimacao_diagnosticos.csv", replace
restore

preserve
    use "$out/estimacao_sx_detalhe.dta", clear
    sort serie disc rel_time cohort
    export delimited using "$out/estimacao_sx_detalhe.csv", replace
restore

*===========================================================================
* 6. Hausman test - S&X aproximado vs DR-DiD
*===========================================================================

di _newline(2) "=== HAUSMAN TEST: S&X aproximado vs DR-DiD ==="
di "(usando variancia conservadora = soma de variancias marginais)"

foreach serie in "9EF" "EM3" {
foreach cocomp in 1 2 {

    if `cocomp' == 1 local disc "LP"
    if `cocomp' == 2 local disc "Mat"

    capture use "$out/sx_`serie'_`disc'.dta", clear
    if _rc {
        di "  `serie' `disc': sem S&X salvo, pulando."
        continue
    }
    quietly count
    if r(N) == 0 {
        di "  `serie' `disc': S&X vazio (nenhuma estimativa valida), pulando."
        continue
    }

    capture estimates restore drdid_plot_`serie'_`disc'
    if _rc {
        di "  `serie' `disc': sem DR-DiD padronizado armazenado."
        continue
    }

    local H_stat      = 0
    local k_used      = 0
    local worst_l     = .
    local worst_diff  = 0

    local lmin = -$K_PRE
    local lmax =  $K_POST

    forvalues l = `lmin'(1)`lmax' {
        if `l' == -1 continue

        qui sum att if l == `l'
        if r(N) == 0 continue
        local att_sx = r(mean)
        qui sum se if l == `l'
        local se_sx = r(mean)
        if missing(`se_sx') | `se_sx' == 0 continue

        if `l' < 0 {
            local absL = abs(`l')
            capture local att_dr = _b[Tm`absL']
            capture local se_dr  = _se[Tm`absL']
        }
        else {
            capture local att_dr = _b[Tp`l']
            capture local se_dr  = _se[Tp`l']
        }
        if _rc continue
        if missing(`att_dr') | missing(`se_dr') continue

        local diff = `att_sx' - `att_dr'
        local se_diff_cons = sqrt(`se_sx'^2 + `se_dr'^2)
        if `se_diff_cons' == 0 continue

        local H_stat = `H_stat' + ((`diff') / `se_diff_cons')^2
        local ++k_used

        if abs(`diff') > abs(`worst_diff') {
            local worst_diff = `diff'
            local worst_l    = `l'
        }
    }

    if `k_used' > 0 {
        local p_val = 1 - chi2(`k_used', `H_stat')
        di _newline "  `serie' `disc':"
        di "    H (conservador) = " %9.3f `H_stat' " | df = " %3.0f `k_used' ///
           " | p = " %7.4f `p_val'
        if !missing(`worst_l') {
            di "    Maior diferenca em l=" %3.0f `worst_l' ///
               " -> diff = " %7.4f `worst_diff'
        }
        if `p_val' < 0.05 {
            di "    -> Rejeita H0 com variancia conservadora:"
            di "       mudanca composicional materialmente impacta o ATT."
        }
        else {
            di "    -> Nao rejeita H0 com variancia conservadora."
            di "       (teste conservador sub-rejeita; bootstrap daria p-valor menor)"
        }
    }

}
}

di _newline "=== Estimacao concluida ==="
di "Resultados: $out/"
di "  - estimacao_eventstudy_long.csv / .dta"
di "  - estimacao_resumo_tese.csv / .dta"
di "  - estimacao_diagnosticos.csv / .dta"
di "  - estimacao_sx_detalhe.csv / .dta"
di "  - $estim_log"
di "Figuras:    $fig/"

log close
capture copy "$estim_log" "$out/02_estimacao_console_latest.log", replace
