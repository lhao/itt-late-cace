cd /Users/lingxinhao/Documents/0data/AID-LEARN/1_jun-sep
capture log close
log using a9_pnas.log, replace

* Oct 1, 2026
* a9_pnas.do	 FA26
* PNAS main text
* IV LATE
* no vce, vce
* conditional

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

* conditional IV LATE w/o vce
ivregress 2sls y $X i.week_id (tcomp = r), first

* conditional IV LATE w vce
ivregress 2sls y $X i.week_id (tcomp = r), first vce(cluster user_id)


**************
* mode correct
drop y
gen y=ycorrect

* conditional LATE w/o vce
biprobit (y = tcomp $X i.week_id) (tcomp = r $X)
margins, dydx(tcomp) atmeans predict(pmarg1) force

* conditional LATE vce
biprobit (y = tcomp $X i.week_id) (tcomp = r $X), vce(cluster user_id)
margins, dydx(tcomp) atmeans predict(pmarg1) force

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

* conditional IV LATE w/o vce
ivregress 2sls y $X i.week_id (tcomp = r), first

* conditional IV LATE w vce
ivregress 2sls y $X i.week_id (tcomp = r), first vce(cluster user_id)


**************
* mode correct
drop y
gen y=ycorrect

* conditional LATE w/o vce
biprobit (y = tcomp $X i.week_id) (tcomp = r $X)
margins, dydx(tcomp) atmeans predict(pmarg1) force

* conditional LATE vce
biprobit (y = tcomp $X i.week_id) (tcomp = r $X), vce(cluster user_id)
margins, dydx(tcomp) atmeans predict(pmarg1) force


log close
