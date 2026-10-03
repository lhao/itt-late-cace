cd /Users/lingxinhao/Documents/0data/AID-LEARN/1_jun-sep
capture log close
log using a9_pnas_si.log, replace

* Oct 1, 2026
* a9_pnas_si.do	 FA26
* PNAS SI
* ITT
* IV LATE
* LC CACE
* no vce, vce(cluster user_id)

* dtable
use stat_ana, clear
drop if stem==0
tab arm 
dtable, by(arm) nosample continuous(yscore ChatGPT_exp nqs, statistics(mean sd)) factor(ycorrect bh subject week_id, statistics(fvprop)) nformat(%5.2f)

tab group
dtable, by(group) nosample continuous(yscore ChatGPT_exp nqs, statistics(mean sd)) factor(ycorrect bh subject week_id, statistics(fvprop)) nformat(%5.2f)

*=============
* Arm2 vs Arm1 
*=============
use stat_ana, clear

tab stem
drop if stem==0
drop if arm==3
tab arm group, nol

* r is assignment comp is complier for treatment group
gen r=0 if arm==1
replace r=1 if arm==2
gen comp=.
replace comp=0 if group==2
replace comp=1 if group==3 
tab r comp, m

* tcomp includes control complied
gen tcomp=comp 
replace tcomp=0 if comp==.

* naming covariates
tab week_id, gen(week)
tab subject, gen(sub)

global X ChatGPT_exp bh sub2-sub5 nqs 

************
* mean score
gen y=yscore

* conditional ITT using OLS
reg y r $X
est store m1a

reg y r $X, vce(cluster user_id)
est store m2a

* conditional IV LATE w/o vce
ivregress 2sls y $X i.week_id (tcomp = r), first
est store m3a

* conditional IV LATE w vce
ivregress 2sls y $X i.week_id (tcomp = r), first vce(cluster user_id)
est store m4a

* conditional CACE ER w/o vce
# delimit ;
gsem 
	(1.C: y <- i.r@0)
	(2.C: y <- i.r) 
	(y <- $X i.week_id)
	(C <- $X)
	(1: comp <- _cons@-15, logit) 
	(2: comp <- _cons@15, logit)
	, lclass(C 2) nolog lcinvariant(coef) 
;
# delimit cr
est store m5a

* conditional CACE ER w vce
# delimit ;
gsem 
	(1.C: y <- i.r@0)
	(2.C: y <- i.r) 
	(y <- $X i.week_id)
	(C <- $X)
	(1: comp <- _cons@-15, logit) 
	(2: comp <- _cons@15, logit)
	, lclass(C 2) nolog lcinvariant(coef) vce(cluster user_id)
;
# delimit cr
est store m6a

**************
* mode correct
drop y
gen y=ycorrect

* conditional ITT using probit
probit y r $X
est store m1b

probit y r $X, vce(cluster user_id)
est store m2b

* conditional LATE w/o vce
biprobit (y = tcomp $X i.week_id) (tcomp = r $X)
est store m3b
margins, dydx(tcomp) atmeans predict(pmarg1) force

* conditional LATE vce
biprobit (y = tcomp $X i.week_id) (tcomp = r $X), vce(cluster user_id)
est store m4b
margins, dydx(tcomp) atmeans predict(pmarg1) force

* conditional CACE ER w/o vce
# delimit ;
gsem 
  (1.C: y <- i.r@0  , family(bernoulli) link(probit)) 
  (2.C: y <- i.r    , family(bernoulli) link(probit)) 
  (y <- $X i.week_id, family(bernoulli) link(probit)) 
	(C <- $X)
	(1:comp <- _cons@-15, logit) 
	(2:comp <- _cons@15, logit)
	, lclass(C 2) nolog lcinvariant(coef)
;
# delimit cr
est store m5b
margins, dydx(r) predict(outcome(y) class(2)) atmeans force

* conditional CACE ER w vce
# delimit ;
gsem 
  (1.C: y <- i.r@0  , family(bernoulli) link(probit)) 
  (2.C: y <- i.r    , family(bernoulli) link(probit)) 
  (y <- $X i.week_id, family(bernoulli) link(probit)) 
	(C <- $X)
	(1:comp <- _cons@-15, logit) 
	(2:comp <- _cons@15, logit)
	, lclass(C 2) nolog lcinvariant(coef) vce(cluster user_id)
;
# delimit cr
est store m6b
margins, dydx(r) predict(outcome(y) class(2)) atmeans force

*********************
* Arm2 vs Arm1 output
* hand-enter marginal effects (6 cols)
* c-ITT, c-LATE c-CACE

* c-ITT, c-LATE
etable, estimates(m2a m4a m2b m4b ) showstars mstat(ll) mstat(rank) mstat(aic) mstat(bic) mstat(N) export(a9_pnas_si_12_coef1.docx, replace)

* c-CACE
etable, estimates(m6a m6b) showstars mstat(ll) mstat(rank) mstat(aic) mstat(bic) mstat(N) export(a9_pnas_si_12_coef2.docx, replace)


*=============
* Arm3 vs Arm1 
*=============
use stat_ana, clear
tab stem
drop if stem==0
drop if arm==2
tab arm group, nol

gen r=0 if arm==1
replace r=1 if arm==3
gen comp=.
replace comp=0 if group==4
replace comp=1 if group==5
tab r comp, m

* tcomp includes control complied
gen tcomp=comp 
replace tcomp=0 if comp==.

* naming covariates
tab week_id, gen(week)
tab subject, gen(sub)

global X ChatGPT_exp bh sub2-sub5 nqs 

************
* mean score
gen y=yscore

* conditional ITT using OLS
reg y r $X
est store n1a

reg y r $X, vce(cluster user_id)
est store n2a

* conditional IV LATE w/o vce
ivregress 2sls y $X i.week_id (tcomp = r), first
est store n3a

* conditional IV LATE w vce
ivregress 2sls y $X i.week_id (tcomp = r), first vce(cluster user_id)
est store n4a

* conditional CACE ER w/o vce
# delimit ;
gsem 
	(1.C: y <- i.r@0)
	(2.C: y <- i.r) 
	(y <- $X i.week_id)
	(C <- $X)
	(1: comp <- _cons@-15, logit) 
	(2: comp <- _cons@15, logit)
	, lclass(C 2) nolog lcinvariant(coef) 
;
# delimit cr
est store n5a

* conditional CACE ER w vce
# delimit ;
gsem 
	(1.C: y <- i.r@0)
	(2.C: y <- i.r) 
	(y <- $X i.week_id)
	(C <- $X)
	(1: comp <- _cons@-15, logit) 
	(2: comp <- _cons@15, logit)
	, lclass(C 2) nolog lcinvariant(coef) vce(cluster user_id)
;
# delimit cr
est store n6a

**************
* mode correct
drop y
gen y=ycorrect

* conditional ITT using probit
probit y r $X
est store n1b

probit y r $X, vce(cluster user_id)
est store n2b

* conditional LATE w/o vce
biprobit (y = tcomp $X i.week_id) (tcomp = r $X)
est store n3b
margins, dydx(tcomp) atmeans predict(pmarg1) force

* conditional LATE vce
biprobit (y = tcomp $X i.week_id) (tcomp = r $X), vce(cluster user_id)
est store n4b
margins, dydx(tcomp) atmeans predict(pmarg1) force

* conditional CACE ER w/o vce
# delimit ;
gsem 
  (1.C: y <- i.r@0  , family(bernoulli) link(probit)) 
  (2.C: y <- i.r    , family(bernoulli) link(probit)) 
  (y <- $X i.week_id, family(bernoulli) link(probit)) 
	(C <- $X)
	(1:comp <- _cons@-15, logit) 
	(2:comp <- _cons@15, logit)
	, lclass(C 2) nolog lcinvariant(coef)
;
# delimit cr
est store n5b
margins, dydx(r) predict(outcome(y) class(2)) atmeans force

* conditional CACE ER w vce
# delimit ;
gsem 
  (1.C: y <- i.r@0  , family(bernoulli) link(probit)) 
  (2.C: y <- i.r    , family(bernoulli) link(probit)) 
  (y <- $X i.week_id, family(bernoulli) link(probit)) 
	(C <- $X)
	(1:comp <- _cons@-15, logit) 
	(2:comp <- _cons@15, logit)
	, lclass(C 2) nolog lcinvariant(coef) vce(cluster user_id)
;
# delimit cr
est store n6b
margins, dydx(r) predict(outcome(y) class(2)) atmeans force

*********************
* Arm3 vs Arm1 output
* hand-enter marginal effects (6 cols)
* c-ITT, c-LATE c-CACE

* c-ITT, c-LATE
etable, estimates(n2a n4a n2b n4b ) showstars mstat(ll) mstat(rank) mstat(aic) mstat(bic) mstat(N) export(a9_pnas_si_13_coef1.docx, replace)

* c-CACE
etable, estimates(n6a n6b) showstars mstat(ll) mstat(rank) mstat(aic) mstat(bic) mstat(N) export(a9_pnas_si_13_coef2.docx, replace)

