/*===========================================================================
  01b_merge_controles.do - Adiciona covariadas ao painel escola-ano
  Andreas Azambuja Barbisan - TCC Insper 2026

  Inputs:
    checkpoints/painel_escola_ano.dta      (gerado por 01_prep.do)
    checkpoints/matricula_escola_ano_serie.dta (gerado por 00_prep_controles.do)
    Agregado Saresp/2011.csv               (proficiencia baseline)

  Output:
    checkpoints/painel_completo.dta
===========================================================================*/

global root   "d:/Andreas/Teste Claude/TCC"
global saresp "$root/Agregado Saresp"
global out    "$root/checkpoints"

use "$out/painel_escola_ano.dta", clear

*---------------------------------------------------------------------------
* 1. MERGE: covariadas de matricula (escola-ano-serie)
*    Microdados entram aqui apenas como estatisticas agregadas de composicao.
*---------------------------------------------------------------------------

local mat_covs ///
    n_alunos pct_fem pct_masc pct_nee idade_media ///
    pct_branca pct_preta pct_parda pct_preta_parda pct_raca_declarada ///
    pct_bolsa_fam pct_rural pct_integral

    merge m:1 codesc year serie using "$out/matricula_escola_ano_serie.dta", ///
    keepusing(`mat_covs') keep(master match) nogen

label var n_alunos           "Matriculas estaduais ativas na escola-serie"
label var pct_fem            "Fracao de alunas na escola-serie"
label var pct_masc           "Fracao de alunos homens na escola-serie"
label var pct_nee            "Fracao de alunos com deficiencia registrada na escola-serie"
label var idade_media        "Idade media dos matriculados na escola-serie"
label var pct_preta_parda    "Fracao preta ou parda na escola-serie"
label var pct_raca_declarada "Fracao com raca/cor declarada na escola-serie"
label var pct_bolsa_fam      "Fracao com Bolsa Familia na escola-serie"

*---------------------------------------------------------------------------
* 1b. MERGE: controles municipais
*---------------------------------------------------------------------------

merge m:1 cd_ibge6 year using "$out/controles_municipais.dta", ///
    keepusing(pib pop pib_pc ln_pib_pc ln_pop prop_jovem_0_19) keep(master match) nogen

label var ln_pib_pc       "Log do PIB per capita municipal"
label var ln_pop          "Log da populacao municipal"
label var prop_jovem_0_19 "Proporcao municipal da populacao de 0 a 19 anos"

*---------------------------------------------------------------------------
* 2. PROFICIENCIA BASELINE (2011) por escola-serie-disciplina
*---------------------------------------------------------------------------

preserve
    import delimited "$saresp/2011.csv", ///
        delimiter(";") encoding("latin1") clear case(lower) varnames(1)

    capture rename medprof MEDPROF
    capture rename MEDPROF medprof
    replace medprof = subinstr(medprof, ",", ".", .)
    destring medprof, replace force

    keep if depadm == 1
    keep if inlist(co_comp, 1, 2)
    drop if missing(medprof)

    gen serie = ""
    replace serie = "3EF"  if strpos(serie_ano, "3")  & strpos(serie_ano, "Ano")
    replace serie = "5EF"  if strpos(serie_ano, "5")  & strpos(serie_ano, "Ano")
    replace serie = "7EF"  if strpos(serie_ano, "7")  & strpos(serie_ano, "Ano")
    replace serie = "9EF"  if strpos(serie_ano, "9")  & strpos(serie_ano, "Ano")
    replace serie = "EM3"  if strpos(serie_ano, "EM") | strpos(serie_ano, "rie")
    drop if serie == ""

    collapse (mean) medprof_base=medprof, by(codesc serie co_comp)
    tempfile base2011
    save `base2011'
restore

merge m:1 codesc serie co_comp using `base2011', keep(master match) nogen
label var medprof_base "Proficiencia media 2011 (baseline)"

*---------------------------------------------------------------------------
* 3. COVARIADAS DERIVADAS E IMPUTACAO CONSERVADORA
*    Mantemos variaveis brutas para diagnostico. As versoes *_imp entram na
*    estimacao para evitar que poucos missings derrubem a amostra inteira.
*---------------------------------------------------------------------------

gen byte miss_n_alunos = missing(n_alunos)
gen ln_alunos_raw = log(n_alunos) if n_alunos > 0
gen double n_alunos_imp = n_alunos
bysort year: egen double ymean_n_alunos = mean(n_alunos)
quietly summarize n_alunos, meanonly
local grand_n_alunos = r(mean)
replace n_alunos_imp = ymean_n_alunos if missing(n_alunos_imp)
replace n_alunos_imp = `grand_n_alunos' if missing(n_alunos_imp) & !missing(`grand_n_alunos')
drop ymean_n_alunos

gen ln_alunos = log(n_alunos_imp) if n_alunos_imp > 0
label var n_alunos_imp "Matriculas estaduais ativas (imputado por media do ano)"
label var miss_n_alunos "Indicador de n_alunos ausente no merge"
label var ln_alunos_raw "Log de matriculas estaduais ativas sem imputacao"
label var ln_alunos "Log de matriculas estaduais ativas"

foreach v in pct_fem pct_nee idade_media pct_preta_parda pct_raca_declarada pct_bolsa_fam {
    gen byte miss_`v' = missing(`v')
    gen double `v'_imp = `v'
    bysort year: egen double ymean_`v' = mean(`v')
    quietly summarize `v', meanonly
    local grand_`v' = r(mean)
    replace `v'_imp = ymean_`v' if missing(`v'_imp)
    replace `v'_imp = `grand_`v'' if missing(`v'_imp) & !missing(`grand_`v'')
    drop ymean_`v'
    label var miss_`v' "Indicador de `v' ausente no merge"
    label var `v'_imp "`v' imputado por media do ano"
}

capture confirm variable zona
if !_rc {
    gen zona_id = zona
}
else {
    encode descricao_zona, gen(zona_id)
}

capture confirm variable tipoesc
if !_rc {
    gen tipo_id = tipoesc
}
else {
    encode descricao_tipo_escola, gen(tipo_id)
}

*---------------------------------------------------------------------------
* 4. SAVE
*---------------------------------------------------------------------------

sort codesc year serie co_comp
save "$out/painel_completo.dta", replace

*---------------------------------------------------------------------------
* 5. DIAGNOSTICO
*---------------------------------------------------------------------------

di ""
di "=== PAINEL COMPLETO ==="
di "Obs: `=_N'"
di ""
di "Cobertura de covariadas brutas:"
foreach v in pct_fem pct_nee idade_media pct_preta_parda pct_raca_declarada pct_bolsa_fam n_alunos medprof_base {
    qui count if missing(`v')
    di "  `v' missing: `r(N)' (`= round(r(N)/_N*100, .1)'%)"
}

di ""
di "Covariadas de composicao usadas na estimacao:"
di "  pct_fem_imp pct_nee_imp idade_media_imp pct_preta_parda_imp pct_raca_declarada_imp pct_bolsa_fam_imp ln_alunos"

di ""
di "Pronto. Arquivo salvo em: $out/painel_completo.dta"
di "Proximo passo: rodar 02_estimacao.do"

