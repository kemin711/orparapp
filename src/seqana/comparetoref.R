/**
 The cutoff defined in input
 */

# calculate fraction of alignments > cutoff 
# from the peak region
aln.files <- c("April27/rawlib.info", "May12/rawlib.info",
	"May13/rawlib.info", "May16/rawlib.info",
"May16/turnoff3trimming.info", "June25/rawlib.info", "July9/rawlib.info",
"July12/rawlib.info", "July15/rawlib.info");

avg <- function(alndf) {
	tmp <- sum(alndf$identity*alndf$count);
	tmp2 <- sum(alndf$coverage*alndf$count);
	total <- sum(alndf$count)
	return(c(tmp/total, tmp2/total));
}

countaln <- function(idencut, covcut) {
	pkaln <- rep(NA,9);
	pkaln.cnt <- rep(NA, 9);
	avgiden <- rep(NA,9);
	avgcov <- rep(NA,9);
	for (i in 1:length(aln.files)) {
		df <- read.table(aln.files[i], header=T);
		pr <- peaks[i,];
		dfp <- subset(df, readlen > pr[1] & readlen < pr[2]);
		pcnt <- sum(dfp$count); #num seq in peak
		dfp.good <- subset(dfp, identity>idencut & coverage>covcut);
		dfp.good.cnt <- sum(dfp.good$count); # num good aln in peak
		avgfeat <- avg(dfp)
		print(paste(pcnt, '  ', dfp.good.cnt) );
		print(paste("average identity, coverage: ", avgfeat[1], ', ', avgfeat[2]));
		pkaln[i] <- round(dfp.good.cnt/pcnt, 4);
		pkaln.cnt[i] <- dfp.good.cnt;
		avgiden[i] <- avgfeat[1];
		avgcov[i] <- avgfeat[2];
	}
	return(list(pkalnfrac=pkaln, pkalncnt=pkaln.cnt), avgiden=avgiden, avgcov=avgcov);
}

res85 <- countaln(0.85, 0.85)
res90 <- countaln(0.9, 0.9)
res8885 <- countaln(0.88, 0.85)
res9186 <- countaln(0.91, 0.86)


print(pkaln)
print(pkaln.cnt)


countGoodf <- function(filepath, idencut, covcut) {
	df <- read.table(filepath, header=T);
	totalcnt <- sum(df$count);
	dfg <- subset(df, identity > frac & coverage > frac);
	goodcnt <- sum(dfg$count);
	x <- list(good=goodcnt,total=totalcnt)
	return (x)
}

input.files <- ("April27/rawlib.info", "May13/rawlib.info", "May16/default.info",
"May16/rawlib.info", "May12/rawlib.info", "June25/rawlib.info", "July9/rawlib.info",
"July12/rawlib.info", "July15/rawlib.info");

countGood("April27/rawlib.info");
countGood("May13/rawlib.info");
countGood("May16/default.info");
countGood("May16/rawlib.info");
countGood("May12/rawlib.info");
countGood("June25/rawlib.info");
countGood("July9/rawlib.info");
countGood("July12/rawlib.info");
countGood("July15/rawlib.info");

countGoodf("June25/rawlib.info", 0.89);

df1 <- read.table("June25/rawlib.info", header=T);
df2 <- read.table("July9/rawlib.info", header=T);

df1g <- subset(df1, identity>0.89 & coverage>0.89);
df2g <- subset(df2, identity>0.89 & coverage>0.89);


goodLength1 <- subset(df1g, readlen > 340 & readlen<400)
sum(goodLength1$count);
sum(df1g$count)

goodLength2 <- subset(df2g, readlen > 340 & readlen<400)
sum(goodLength2$count);
sum(df2g$count)


/*
 length distribution
*/

lcount1 <- read.table("June25/length.dis", header=T);
lcount2 <- read.table("July9/length.dis", header=T);
xlimit <- c(50, 500)

par(mfrow=c(2,1))
plot(lcount1$length, lcount1$count, xlim=xlimit);
plot(lcount2$length, lcount2$count, xlim=xlimit);

g1 <- subset(lcount1, length>340 & length < 400)
g2 <- subset(lcount2, length>340 & length < 400)

sum(g1$count)/sum(lcount1$count)
sum(g2$count)/sum(lcount2$count)


-- this will count the good and total alignments
countGood <- function(filepath) {
	df <- read.table(filepath, header=T);
	totalcnt <- sum(df$count);
	dfg <- subset(df, identity > 0.85 & coverage > 0.85);
	goodcnt <- sum(dfg$count);
	x <- list(good=goodcnt,total=totalcnt)
	return (x)
}

countGood("April27/rawlib.info");
countGood("May13/rawlib.info");
countGood("May16/default.info");
countGood("May16/rawlib.info");
countGood("June25/rawlib10percent.info");
countGood("May12/rawlib.info");
countGood("June25/rawlib.info");
countGood("July9/rawlib.info");
countGood("July12/rawlib.info");
countGood("July15/rawlib.info");

/* result
June25
	short   1043590
	artifacts  1080 
	rest    6168875 
   Alignment:
	good   175791
	total 3421508

July9
	short 1390357
	artifacts 285 
	rest  5236819 
   Alignment:
	good     2791
	total 2073541

July12
	short 75141 short sequences
	artifacts 171 
	rest 5466050 
   Alignment:
	good        9
	total 4136052

July15
	5015975 short sequences
	569 artifacts
	2483336 rest
   Alignment:
	good    70239
	total 1342857
*/
-- on Kraken /scratch2/zhouke/virology
countGood("April27/rawlib.info");
countGood("May13/rawlib.info");
countGood("May16/default.info");
countGood("May16/rawlib.info");

countGood("June25/rawlib10percent.info");

countGoodf <- function(filepath, frac) {
	df <- read.table(filepath, header=T);
	totalcnt <- sum(df$count);
	dfg <- subset(df, identity > frac & coverage > frac);
	goodcnt <- sum(dfg$count);
	x <- list(good=goodcnt,total=totalcnt)
	return (x)
}

countGoodf("June25/test.info", 0.89);

/*
result
April27
	128317 short sequences
	349 artifacts
	6206600 rest
   alignment:
	good  2586861
	total 5513848

May13
	533984 short sequences
	530 artifacts
	6169737 rest
	Alignment:
	good  1074661
	total 2995987

May16
	Default
	 224655 short sequences
	    710 artifacts
	6540041 rest

	Alignment:
	good  1,205,069
	total 4,041,108

	Rawlib
	224655 short sequences
	710 artifacts
	6540041 rest
	Alignment
	good  1205069
	total 4041108

	turn off 3 trimming
	69562 short sequences
	1480 artifacts
	6791879 rest
	Alignment:
	still running, quality much worse than the above 2,

*/



sn131 <- read.table("viroSN1_31.info", header=T);
sn133 <- read.table("SN1_33_system92.info", header=T);
sn132 <- read.table("virobarcode_91_070.info", header=T);

xrange <- c(2800,3250)
par(mfrow=c(2,2))
plot(sn131$score, sn131$count, xlim=xrange);
plot(sn132$score, sn132$count, xlim=xrange, ylim=c(0,1500));
plot(sn133$score, sn133$count, xlim=xrange, ylim=c(0,9000));

# full view of the same height
xrange <- c(0,3250)
yrange <- c(0, 100000)
par(mfrow=c(2,2))
plot(sn131$score, sn131$count, xlim=xrange, ylim=yrange);
plot(sn132$score, sn132$count, xlim=xrange, ylim=yrange);
plot(sn133$score, sn133$count, xlim=xrange, ylim=yrange);

# full view but automatic height
xrange <- c(0,3250)
yrange <- c(0, 100000)
par(mfrow=c(2,2))
plot(sn131$score, sn131$count, xlim=xrange);
plot(sn132$score, sn132$count, xlim=xrange);
plot(sn133$score, sn133$count, xlim=xrange);

--identity
par(mfrow=c(2,2))
plot(sn131$identity, sn131$count);
plot(sn132$identity, sn132$count);
plot(sn133$identity, sn133$count);

par(mfrow=c(2,2))
plot(sn131$coverage, sn131$count);
plot(sn132$coverage, sn132$count);
plot(sn133$coverage, sn133$count);


idcut <- 0.88
covcut <- 0.88
--identity
par(mfrow=c(2,2))
for (i in 1:3) {
	sub <- subset(ds[[i]], identity>idcut & coverage>covcut)
	plot(sub$identity, sub$count);
}

par(mfrow=c(2,2))
for (i in 1:3) {
	sub <- subset(ds[[i]], identity>idcut & coverage>covcut)
	plot(sub$coverage, sub$count);
	print(sum(sub$count)/sum(ds[[i]]))
}


par(mfrow=c(2,2))
plot(sn131$coverage, sn131$count);
plot(sn132$coverage, sn132$count);
plot(sn133$coverage, sn133$count);

numaln <- c(nrow(sn131), nrow(sn132), nrow(sn133));

ds <- list(sn131, sn132, sn133);

result <- numeric(0);

for (d in ds) {
   nr <- nrow(subset(d, identity > 0.85 & coverage > 0.8));
	frac <- nr/nrow(d)
	frac
}



nrow(subset(sn131, identity>idcut & coverage>covcut))/nrow(sn131)
nrow(subset(sn132, identity>idcut & coverage>covcut))/nrow(sn132)
nrow(subset(sn133, identity>idcut & coverage>covcut))/nrow(sn133)

# nrow identity>0.85 & coverage>0.8
sn131 5517
sn132 125
sn133 6417

meaniden <- function(df) {
	avg <- sum((df$identity*df$count))/sum(df$count);
	avg;
}
meancov <- function(df) {
	avg <- sum((df$coverage*df$count))/sum(df$count);
	avg;
}

# got overflow
meanscore <- function(df) {
	avg <- sum((df$score*df$count))/sum(df$count);
	avg;
}

for (d in ds) {
   s <- subset(d, score > 2600);
	frac <- sum(s$count)/sum(d);
	print(frac);
	print(paste(' identity: ', meaniden(d)));
	print(paste(' coverage: ', meancov(d)));
}




