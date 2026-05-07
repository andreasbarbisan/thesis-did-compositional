/*===========================================================================
  01_prep.do — Preparação do painel escola-ano para análise DiD PEI
  Andreas Azambuja Barbisan — TCC Insper 2026

  Inputs:
    Agregado Saresp/XXXX.csv  (2011–2018)
    TCC v2 Micro 4/bruto/escolas_pei_first_treat.dta

  Output:
    checkpoints/painel_escola_ano.dta   — painel limpo, pronto para estimação
===========================================================================*/

clear all
set more off

global root   "d:/Andreas/Teste Claude/TCC"
global saresp "$root/Agregado Saresp"
global pei    "$root/TCC v2 Micro 4/bruto/escolas_pei_first_treat.dta"
global ender  "$root/checkpoints/endereco_escolas.dta"
global out    "$root/checkpoints"

capture mkdir "$out"

*---------------------------------------------------------------------------
* 1. IMPORTAR E EMPILHAR OS CSVs AGREGADOS (2011–2018)
*---------------------------------------------------------------------------

tempfile empilhado
local primeiro = 1

foreach ano in 2011 2012 2013 2014 2015 2016 2017 2018 {

    import delimited "$saresp/`ano'.csv", ///
        delimiter(";") encoding("latin1") clear case(lower) varnames(1)

    * Padronizar colunas-chave (casing varia entre anos)
    capture rename codesc    codesc
    capture rename serie_ano serie_ano
    capture rename co_comp   co_comp
    capture rename ds_comp   ds_comp
    capture rename depadm    depadm
    capture rename depbol    depbol

    * MEDPROF/medprof vem como string com vírgula decimal
    capture rename medprof MEDPROF
    capture rename MEDPROF medprof
    replace medprof = subinstr(medprof, ",", ".", .)
    destring medprof, replace force

    keep codesc serie_ano co_comp ds_comp medprof depadm depbol
    gen year = `ano'

    if `primeiro' {
        save `empilhado', replace
        local primeiro = 0
    }
    else {
        append using `empilhado'
        save `empilhado', replace
    }

    di "  -> `ano': OK (`=_N' obs acumuladas)"
}

*---------------------------------------------------------------------------
* 2. FILTROS E HARMONIZAÇÃO
*---------------------------------------------------------------------------

* Manter só rede estadual (PEI é programa estadual)
keep if depadm == 1

* Manter só LP (co_comp==1) e Matemática (co_comp==2)
keep if inlist(co_comp, 1, 2)
label define disc 1 "LP" 2 "Matematica"
label values co_comp disc

* Harmonizar série sem depender de acentos (encoding pode corromper)
gen serie = ""
replace serie = "3EF"  if strpos(serie_ano, "3")  & strpos(serie_ano, "Ano")
replace serie = "5EF"  if strpos(serie_ano, "5")  & strpos(serie_ano, "Ano")
replace serie = "7EF"  if strpos(serie_ano, "7")  & strpos(serie_ano, "Ano")
replace serie = "9EF"  if strpos(serie_ano, "9")  & strpos(serie_ano, "Ano")
replace serie = "EM3"  if strpos(serie_ano, "EM") | strpos(serie_ano, "rie")

drop if serie == ""
drop serie_ano ds_comp depbol
drop if missing(medprof)

* Colapsar para único valor por escola-ano-série-disciplina
* (elimina duplicatas por turno/período dentro da escola)
collapse (mean) medprof, by(codesc year serie co_comp depadm)

*---------------------------------------------------------------------------
* 3. MERGE COM CADASTRO DE ESCOLAS E BASE PEI
*---------------------------------------------------------------------------

merge m:1 codesc using "$ender", ///
    keepusing(cd_ibge7 cd_ibge6 de mun distr latitude longitude zona tipoesc nomesc) ///
    keep(master match) nogen

merge m:1 codesc using "$pei", ///
    keepusing(first_treat dre descricao_tipo_escola descricao_zona carga_horaria) ///
    keep(master match)

replace first_treat = 0 if _merge == 1
drop _merge

*---------------------------------------------------------------------------
* 4. VARIÁVEIS DiD
*---------------------------------------------------------------------------

gen treated   = (first_treat > 0) if !missing(first_treat)

gen post      = (year >= first_treat) if treated == 1
replace post  = 0 if treated == 0

gen rel_time  = year - first_treat if treated == 1
replace rel_time = . if treated == 0

* gvar: grupo para csdid (= ano de adesão; 0 para never-treated)
gen gvar = first_treat
replace gvar = 0 if treated == 0

*---------------------------------------------------------------------------
* 5. LABELS E SAVE
*---------------------------------------------------------------------------

label var codesc      "Codigo da escola"
label var year        "Ano"
label var serie       "Serie avaliada"
label var co_comp     "Disciplina (1=LP 2=Mat)"
label var medprof     "Proficiencia media SARESP"
label var first_treat "Ano de adesao ao PEI (0=nunca)"
label var treated     "Ever-treated (PEI)"
label var post        "Pos-tratamento"
label var rel_time    "Tempo relativo ao tratamento"
label var gvar        "Grupo (para csdid)"
label var mun         "Municipio"
label var de          "Diretoria de ensino"
label var dre         "Diretoria regional de ensino"
label var distr       "Distrito"
label var cd_ibge7    "Codigo IBGE municipal com 7 digitos"
label var cd_ibge6    "Codigo IBGE municipal com 6 digitos"

order codesc year serie co_comp medprof treated gvar first_treat post rel_time ///
      cd_ibge7 cd_ibge6 mun de dre distr latitude longitude ///
      descricao_tipo_escola descricao_zona carga_horaria depadm

sort codesc year serie co_comp

save "$out/painel_escola_ano.dta", replace

*---------------------------------------------------------------------------
* 6. DIAGNÓSTICO RÁPIDO
*---------------------------------------------------------------------------

di ""
di "=== PAINEL FINAL ==="
tab year
tab serie
tab co_comp
tab gvar if treated == 1

di ""
di "Pronto. Arquivo salvo em: $out/painel_escola_ano.dta"
di "Proximo passo: rodar 02_estimacao.do"
