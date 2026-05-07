/*===========================================================================
  02b_sx_nivel_b.do - Implementacao paramétrica do estimador Sant'Anna & Xu
  (2026), nivel B. Substitui o bloco SX2026 original de 02_estimacao.do
  (linhas 351-545). Tudo que nao esta neste arquivo permanece igual.

  Diferencas chave em relacao ao bloco original:
    1. PS multinomial 4-cat (mlogit), nao logit binario.
    2. Tres OR regressions separadas para m_{1,0}, m_{0,1}, m_{0,0}.
    3. Monta tau_dr via equacao (2.6) do artigo, em vez de chamar 'drdid'.
    4. Covariadas fixadas em g-1 (baseline da coorte), evitando bad controls.
    5. Variancia agregada via funcao de influencia empirica, agrupada por
       escola (cluster).
    6. Fallback para logit binario + S&Z se mlogit nao convergir.

  Notacao (equacao 2.6 do artigo):
    tau_dr = E[ w_{1,1}*(Y - (m_{1,0}+m_{0,1}-m_{0,0}))
                + w_{1,0}*(Y - m_{1,0})*(-1)
                - w_{0,1}*(Y - m_{0,1})
                + w_{0,0}*(Y - m_{0,0}) ]

  onde w_{d,t}(D,T,X) = I_{d,t} * p(1,1,X) / p(d,t,X), normalizados Hajek.

  Este bloco assume que existe uma variavel 'cohort_base_*' criada por escola
  com os valores das covariadas em g-1. O codigo cria essa variavel dentro
  do loop para cada coorte ativa.
===========================================================================*/

*-------------------------------------------------------------------
* 4. S&X aproximado (nivel B): mlogit 4-cat + 3 OR separadas
*-------------------------------------------------------------------
di _newline "--- 4. S&X aproximado (nivel B) ---"
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

* Covariadas para o PS multinomial e para as OR.
* Usamos as mesmas sete covariadas composicionais, mas fixadas no valor
* da escola em g-1 (baseline da coorte). Isso evita bad controls.
local cov_sx_list "pct_fem_imp pct_nee_imp idade_media_imp pct_preta_parda_imp pct_raca_declarada_imp pct_bolsa_fam_imp ln_alunos"

local row = 1
forvalues l = `lmin'(1)`lmax' {
    matrix sx_l[`row', 1] = `l'

    if `l' == -1 {
        di "      l=" %3.0f `l' " -> baseline universal normalizado para zero"
        foreach g of local cohorts {
            local t = `g' + `l'
            if `t' < 2011 | `t' > 2018 continue
            post `post_sx' ("`serie'") ("`disc'") (`l') (`g') (`t') ///
                (.) (.) (.) (.) (.) (.) (0) (.) (.) ("baseline_km1")
        }
        local row = `row' + 1
        continue
    }

    local num_w_att = 0    // acumula numerador: soma_g [ att_g * n_11_g ]
    local den_w     = 0    // acumula denominador: soma_g [ n_11_g ]
    local var_agg   = 0    // acumula variancia agregada
    local sx_pairs_ok = 0

    foreach g of local cohorts {
        local t = `g' + `l'
        if `t' < 2011 | `t' > 2018 continue

        qui count if (gvar == `g' | gvar == 0) ///
            & inlist(year, `g' - 1, `t') & !missing(medprof_std)
        local n_pair = r(N)
        if `n_pair' < 30 {
            local ++sx_lowobs_total
            post `post_sx' ("`serie'") ("`disc'") (`l') (`g') (`t') ///
                (`n_pair') (.) (.) (.) (.) (.) (0) (.) (.) ("too_few_obs")
            continue
        }

        local att_g = .
        local se_g  = .
        local n_11  = 0
        local n_10  = .
        local n_01  = .
        local n_00  = .
        local converged = 0
        local status "init"

        capture frame drop gt_frame
        frame copy default gt_frame
        frame change gt_frame

            keep if (gvar == `g' | gvar == 0) & inlist(year, `g' - 1, `t')

            * --- BASELINE COVARIATES: fixa covariadas em g-1 ---
            * Para cada escola, pega o valor das covariadas composicionais
            * no ano g-1 e replica para o ano t. Isso evita usar covariadas
            * pos-tratamento (bad controls).
            foreach v of local cov_sx_list {
                gen double `v'_b = `v' if year == `g' - 1
                bysort codesc (year): egen double `v'_base = max(`v'_b)
                drop `v'_b
            }

            * Droppar escolas sem baseline em g-1 (escolas que so aparecem no pos)
            egen byte miss_base = rowmiss(pct_fem_imp_base pct_nee_imp_base ///
                idade_media_imp_base pct_preta_parda_imp_base ///
                pct_raca_declarada_imp_base pct_bolsa_fam_imp_base ln_alunos_base)
            keep if miss_base == 0
            drop miss_base

            * Variaveis do setup 2x2
            gen byte D_sx    = (gvar == `g')
            gen byte period  = (year == `t')
            gen byte cat_11  = (D_sx == 1) & (period == 1)
            gen byte cat_10  = (D_sx == 1) & (period == 0)
            gen byte cat_01  = (D_sx == 0) & (period == 1)
            gen byte cat_00  = (D_sx == 0) & (period == 0)

            * Variavel categorica para mlogit: 0=(0,0), 1=(0,1), 2=(1,0), 3=(1,1)
            gen byte dtcat = .
            replace dtcat = 0 if cat_00 == 1
            replace dtcat = 1 if cat_01 == 1
            replace dtcat = 2 if cat_10 == 1
            replace dtcat = 3 if cat_11 == 1

            qui count if !missing(medprof_std, dtcat)
            local n_pair = r(N)
            qui count if dtcat == 3
            local n_11 = r(N)
            qui count if dtcat == 2
            local n_10 = r(N)
            qui count if dtcat == 1
            local n_01 = r(N)
            qui count if dtcat == 0
            local n_00 = r(N)

            * Precisamos das quatro celulas nao-vazias com minimo razoavel
            local n_min = min(`n_11', `n_10', `n_01', `n_00')
            if `n_11' < 5 | `n_10' < 5 | `n_01' < 5 | `n_00' < 5 {
                local status "not_2x2_min5"
                local att_g = .
                local se_g  = .
            }
            else {

                * ---------------------------------------------------------
                * (a) PS MULTINOMIAL: mlogit com 4 categorias, base=(0,0)
                * ---------------------------------------------------------
                local cov_base_list ""
                foreach v of local cov_sx_list {
                    local cov_base_list "`cov_base_list' `v'_base"
                }

                capture mlogit dtcat `cov_base_list', base(0) iterate(100)
                local rc_ml = _rc

                if !`rc_ml' & e(converged) == 1 {
                    * Predict probabilities para cada categoria
                    capture predict double p00_hat, outcome(0)
                    capture predict double p01_hat, outcome(1)
                    capture predict double p10_hat, outcome(2)
                    capture predict double p11_hat, outcome(3)

                    * Checar suporte (evitar probabilidades degeneradas)
                    qui summarize p11_hat, meanonly
                    local p11_min = r(min)
                    qui summarize p00_hat, meanonly
                    local p00_min = r(min)
                    qui summarize p01_hat, meanonly
                    local p01_min = r(min)
                    qui summarize p10_hat, meanonly
                    local p10_min = r(min)

                    local ps_ok = (`p11_min' > 0.005 & `p10_min' > 0.005 ///
                                 & `p01_min' > 0.005 & `p00_min' > 0.005)

                    if `ps_ok' {
                        * ---------------------------------------------------------
                        * (b) OR REGRESSIONS: tres regressoes separadas
                        *     m_{1,0}, m_{0,1}, m_{0,0}. Nao precisamos m_{1,1}.
                        * ---------------------------------------------------------
                        capture regress medprof_std `cov_base_list' if cat_10 == 1
                        local rc_or10 = _rc
                        capture predict double m10_hat, xb

                        capture regress medprof_std `cov_base_list' if cat_01 == 1
                        local rc_or01 = _rc
                        capture predict double m01_hat, xb

                        capture regress medprof_std `cov_base_list' if cat_00 == 1
                        local rc_or00 = _rc
                        capture predict double m00_hat, xb

                        if !`rc_or10' & !`rc_or01' & !`rc_or00' {
                            * ---------------------------------------------------------
                            * (c) PESOS HAJEK equacao (2.4) do artigo
                            * ---------------------------------------------------------
                            * Numeradores brutos (antes de normalizacao)
                            gen double w11_raw = cat_11 / 1    // D*T direto
                            gen double w10_raw = cat_10 * p11_hat / p10_hat
                            gen double w01_raw = cat_01 * p11_hat / p01_hat
                            gen double w00_raw = cat_00 * p11_hat / p00_hat

                            * Normalizacao Hajek
                            qui summarize w11_raw, meanonly
                            local E_w11 = r(sum) / r(N)
                            qui summarize w10_raw, meanonly
                            local E_w10 = r(sum) / r(N)
                            qui summarize w01_raw, meanonly
                            local E_w01 = r(sum) / r(N)
                            qui summarize w00_raw, meanonly
                            local E_w00 = r(sum) / r(N)

                            gen double w11 = w11_raw / `E_w11' if `E_w11' > 0
                            gen double w10 = w10_raw / `E_w10' if `E_w10' > 0
                            gen double w01 = w01_raw / `E_w01' if `E_w01' > 0
                            gen double w00 = w00_raw / `E_w00' if `E_w00' > 0

                            * ---------------------------------------------------------
                            * (d) ESTIMAND tau_dr via equacao (2.6)
                            *     tau_dr = E[w11*(Y - m10 - m01 + m00)]
                            *            - E[w10*(Y - m10)]
                            *            + E[w01*(Y - m01)]  <- sinal (-1)^(0+1)=-1
                            *                                  mas w01 entra subtraindo
                            *            + E[w00*(Y - m00)]  <- sinal (-1)^(0+0)=+1
                            *
                            * Nota sobre sinais (eq 2.6):
                            *   sum_(d,t) in - of (-1)^(d+t) * w_{d,t}*(Y - m_{d,t})
                            *   (1,0): sinal (-1)^1 = -1
                            *   (0,1): sinal (-1)^1 = -1
                            *   (0,0): sinal (-1)^0 = +1
                            * ---------------------------------------------------------

                            gen double psi_i = w11 * (medprof_std - m10_hat - m01_hat + m00_hat) ///
                                             - w10 * (medprof_std - m10_hat) ///
                                             - w01 * (medprof_std - m01_hat) ///
                                             + w00 * (medprof_std - m00_hat)

                            qui summarize psi_i, meanonly
                            local att_g = r(mean)

                            * ---------------------------------------------------------
                            * (e) VARIANCIA via funcao de influencia empirica
                            *     Clusterizada em codesc
                            * ---------------------------------------------------------
                            * Centrar: psi_i - tau_hat
                            gen double psi_c = psi_i - `att_g'
                            * Somar dentro de cada escola (cluster)
                            bysort codesc: egen double psi_cluster = total(psi_c)
                            bysort codesc: gen byte school_tag2 = (_n == 1)
                            qui count if school_tag2
                            local n_clust_g = r(N)
                            qui summarize psi_cluster if school_tag2
                            * var = sum_school (sum_t psi_c)^2 / n^2
                            gen double psi_cluster_sq = psi_cluster^2 if school_tag2
                            qui summarize psi_cluster_sq, meanonly
                            local sum_sq = r(sum)
                            qui count
                            local n_total = r(N)
                            local var_g = `sum_sq' / (`n_total'^2)
                            local se_g  = sqrt(`var_g')
                            local converged = 1
                            local status "ok_mlogit"
                        }
                        else {
                            local status "or_fail"
                            local att_g = .
                            local se_g  = .
                        }
                    }
                    else {
                        local status "ps_lowoverlap"
                        local att_g = .
                        local se_g  = .
                    }
                }
                else {
                    * Fallback: mlogit nao convergiu, cai para drdid binario (S&Z)
                    local ++sx_mlogit_fail
                    capture quietly drdid medprof_std `cov_base_list', ///
                        ivar(codesc) time(period) treatment(D_sx) ///
                        dripw cluster(codesc)
                    local rc_bk = _rc
                    if !`rc_bk' {
                        local att_g = e(att1)
                        local se_g  = sqrt(e(attvar1))
                        local converged = 1
                        local status "fallback_sz"
                    }
                    else {
                        local status "mlogit_and_fallback_fail"
                        local att_g = .
                        local se_g  = .
                    }
                }
            }

        frame change default
        frame drop gt_frame

        * Agregacao tipo C&S: peso = n_11 (tratados no pos da coorte)
        if !missing(`att_g') & `n_11' > 0 {
            local num_w_att = `num_w_att' + `att_g' * `n_11'
            local den_w     = `den_w'     + `n_11'
            * Variancia agregada assumindo independencia entre coortes
            * (aproximacao; S&X recomendaria bootstrap sobre IF conjunta,
            * mas isso ja eh uma melhoria sobre o codigo original)
            local var_agg   = `var_agg'   + (`n_11')^2 * (`se_g')^2
            local ++sx_pairs_ok
            local ++sx_pairs_ok_total
        }
        else {
            local ++sx_rc_fail_total
        }

        post `post_sx' ("`serie'") ("`disc'") (`l') (`g') (`t') ///
            (`n_pair') (`n_11') (`n_11') (`n_10') (`n_01') (`n_00') ///
            (`converged') (`att_g') (`se_g') ("`status'")
    }

    if `den_w' > 0 {
        matrix sx_att[`row', 1] = `num_w_att' / `den_w'
        matrix sx_se[`row', 1]  = sqrt(`var_agg') / `den_w'
        di "      l=" %3.0f `l' " -> ATT=" %9.3f (`num_w_att' / `den_w') ///
           " | SE=" %9.3f (sqrt(`var_agg') / `den_w') ///
           " | coortes validas=" %3.0f `sx_pairs_ok'
    }
    else {
        di "      l=" %3.0f `l' " -> sem estimativa valida"
    }
    local row = `row' + 1
}

di "  [SX-aproximado] mlogit falhou em `sx_mlogit_fail' pares coorte-tempo"
di "  [SX-aproximado] baixo N em `sx_lowobs_total' pares"
di "  [SX-aproximado] total de ATT(g,t) validos: `sx_pairs_ok_total'"

*-------------------------------------------------------------------
* Salvar resultados S&X aproximado no mesmo formato do bloco original
*-------------------------------------------------------------------

capture frame drop sx_frame
frame create sx_frame
frame sx_frame {
    clear
    svmat sx_l,   names(l)
    svmat sx_att, names(att)
    svmat sx_se,  names(se)
    rename (l1 att1 se1) (l att se)
    gen ci_low  = att - 1.96 * se
    gen ci_high = att + 1.96 * se
    gen serie_disc = "`serie'_`disc'"
    drop if missing(att)
    save "$out/sx_`serie'_`disc'.dta", replace
}
frame drop sx_frame

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
        local att_now = el(sx_att, `sx_row', 1)
        local se_now  = el(sx_se,  `sx_row', 1)
        if !missing(`att_now') {
            matrix `b_sx'[1, `j'] = `att_now'
            if !missing(`se_now') {
                matrix `V_sx'[`j', `j'] = (`se_now')^2
            }
        }
    }
    local ++j
}

build_sx_est, bmat(`b_sx') vmat(`V_sx') ///
    store(sx_plot_`serie'_`disc') nclust(`n_schools') eventcells(`sx_pairs_ok_total')

post `post_diag' ("`serie'") ("`disc'") ("SXapr") (.) (`n_schools') ///
    (`treated_obs') (`control_obs') (`treated_schools') (`n_cohorts') (`n_schools') (`sx_pairs_ok_total') ///
    ("S&X paramétrico: mlogit 4-cat + 3 OR; cov baseline; mlogit_fail=`sx_mlogit_fail'")
