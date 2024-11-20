# Credits to Henry for this bit of code 

pd62341_mtr=56
pd62341_dep=199

pd63383_mtr=104
pd63383_dep=174


pd62341_mtr/pd62341_dep

qbinom(p=c(0.05,0.95), size=pd62341_dep, prob=(pd62341_mtr/pd62341_dep))[1]/pd62341_dep

qbinom(p=c(0.05,0.95), size=pd63383_dep, prob=(pd63383_mtr/pd63383_dep))/pd63383_dep


?pbinom
