-- use one data frame instead of collecting files from different directoires
-- data source: /home/zhouke/work/nextgen/shrna

-- The number of hits per libarary sequence 
expname <- c( "treat_orig  ", "treat_a118  ", "treat_a228  ", "untreat_orig",
"untreat_a118", "untreat_a228", "first_228   ", "first_T5    ", "first_xh4   ");

nh <- read.table("numhitsdist_postprocess.tab", header=T)

par(mfrow=c(3,3), mar=c(5,4,4,3))
xlimit <- c(-200, 900)
for ( i in 1:9 ) {
	h <- subset(nh, experiment_id==i);
	h1 <- subset(h, h$numhits >= 0);
	h2 <- subset(h, h$numhits < 0);
	low <- min(h$numhits)
	high <- max(h$numhits)
	top <- max(log(h$count))
	top.positive <- max(log(h1$count))
	plot(h2$numhits, log(h2$count), xlab="numhits", ylab="log (count)", xlim=xlimit, col=2)
	points(h1$numhits, log(h1$count), col=3)
	text(mean(xlimit), top*0.9, expname[i])
	text(xlimit[2]*0.7, top*0.8, paste(low, high, sep=" - "))
	abline(v=0)
	abline(h=top.positive)
}

-- plot four representative plots
par(mfrow=c(2,2), mar=c(5,4,4,3))
xlimit <- c(-200, 600)
mvdown <- 0.9
for ( i in c(3, 4, 7, 8) ) {
	h <- subset(nh, experiment_id==i);
	h1 <- subset(h, h$numhits >= 0);
	h2 <- subset(h, h$numhits < 0);
	low <- min(h$numhits)
	high <- max(h$numhits)
	top <- max(log(h$count))
	top.positive <- max(log(h1$count))
	plot(h2$numhits, log(h2$count), xlab="numhits", ylab="log (count)", xlim=xlimit, col=2)
	points(h1$numhits, log(h1$count), col=3)
	if (i == 8) {
		mvdown <- 0.8;
	}
	else {
		mvdown <- 0.9;
	}
	text(mean(xlimit), top*mvdown, expname[i])
	text(xlimit[2]*0.7, top*mvdown*0.9, paste(low, high, sep=" - "))
	abline(v=0)
	abline(h=top.positive)
}

len <- read.table("Length_all.tab", header=T)
factor(len$partial, labels=c("false", "true"))
/*

+----+--------------+
| id | short_name   |
+----+--------------+
|  1 | treat_orig   |
|  2 | treat_a118   |
|  3 | treat_a228   |
|  4 | untreat_orig |
|  5 | untreat_a118 |
|  6 | untreat_a228 |
|  7 | first_228    |
|  8 | first_T5     |
|  9 | first_xh4    |
+----+--------------+

*/


par(mfrow=c(3,3), mar=c(5,4,4,3))
-- partial sequences
xlimit <- c(0, 100)
for ( i in 1:9 ) {
	cl <- subset(len, experiment_id==i & partial==1);
	plot(cl$length, cl$count, type="l", xlim=xlimit, xlab="Length",
		ylab="count");
	mx <- max(cl$count)
	text(xlimit[2]*0.7, mx*0.9, expname[i])
}

-- complete sequences
xlimit <- c(43, 53)
for ( i in 1:9 ) {
	cl <- subset(len, experiment_id==i & partial==0);
	plot(cl$length, cl$count, type="o", xlim=xlimit, xlab="Length",
		ylab="count");
	mx <- max(cl$count)
	text(mean(xlimit)*1.05, mx*0.9, expname[i])
	abline(v=48, col=4, lty=2)
}

-- identity 3D plot
install.packages("scatterplot3d", dependencies = TRUE)

library(scatterplot3d)

/* example command
scatterplot3d(x = data$x, y = data$y, z = data$z)
*/

lri <- read.table("leftright_count.tab", header=T)

jpeg('lridentity_full.jpg')

pdf('lridentity_full.pdf')

par(mfrow=c(3,3), mar=c(5,4,4,3))
for ( i in 1:9 ) {
	iden <- subset(lri, experiment_id==i);
	scatterplot3d(x=iden$left_identity, y=iden$right_identity, z=iden$count,
		xlab="Left", ylab="Right", zlab="Count");
}

dev.off()


pdf('lridentity_left.pdf')
par(mfrow=c(3,3), mar=c(5,4,4,3))
-- ignoring the right 
for ( i in 1:9 ) {
	iden <- subset(lri, experiment_id==i);
	scatterplot3d(x=iden$left_identity, y=iden$right_identity, z=iden$count,
		xlim=c(0.8, 1),
		xlab="Left", ylab="Right", zlab="Count");
}
dev.off()

-- Only look at both left and right matches
pdf('lridentity_corner.pdf')
par(mfrow=c(3,3), mar=c(5,4,4,3))
for ( i in 1:9 ) {
	iden <- subset(lri, experiment_id==i);
	scatterplot3d(x=iden$left_identity, y=iden$right_identity, z=iden$count,
		xlim=c(0.7, 1), ylim=c(0.7, 1),
		xlab="Left", ylab="Right", zlab="Count");
}
dev.off()


-- R should be run above the following directores
dirs=c("shRNA_2_2_8ul_044", "shRNA_T5_protocol_run_047", "user_xh4_reanalyze_045");

filename <- paste(dirs[1], "/identityLeft.tab", sep="");
iden.l <- read.table(filename, header=T);
colnames(iden.l);

-- plotting identity distribution from each run
-- line for left, dots for right
-- identity starts from 72, we can make the starting point lower
savedi <- list();
par(mfrow=c(2,2), mar=c(5,4,4,3))
for (d in dirs) {
	filename <- paste(d, "/identityLeft.tab", sep="");
	rfilename <- paste(d, "/identityRight.tab", sep="");

	iden <- read.table(filename, header=T);
	idenr <- read.table(rfilename, header=T);
	plot(iden[[1]], iden[[2]], type="l", xlab="Identity", ylab="Count", xlim=c(0.85, 1.01));
	points(idenr[[1]], iden[[2]], type="p");
	text(0.85, max(iden$count)*0.8, d, adj=1, pos=4)
	savedi <- c(savedi, iden)
}
plot(savedi[[1]], savedi[[2]]/sum(savedi[[2]]), type='l', xlim=c(0.85,1), ylim=c(0, 0.9), 
	xlab="Identity", ylab="Frequency");
points(savedi[[3]], savedi[[4]]/sum(savedi[[4]]), type='p', pch=2, col="red");
points(savedi[[5]], savedi[[6]]/sum(savedi[[6]]), type='p', pch=3, col="green");


simpleMean <- function(value, count) {
	sum(value*count)/sum(count)
}

simpleStd <- function(value, count) {
	avg <- simpleMean(value, count);
	total <- sum(count);
	std <- sqrt(sum((count/total)*(value - avg)**2))
	c(round(avg, 3), round(std, 3), total)
}

-- length statistics

dirs=c("shRNA_2_2_8ul_044", "shRNA_T5_protocol_run_047", "user_xh4_reanalyze_045");
lengthsave=list();
par(mfrow=c(2,2))
xlimit <- c(0,100)
for (d in dirs) {
	filename <- paste(d, "/insertlencomplete.tab", sep="");
	pfilename <- paste(d, "/insertlenpartial.tab", sep="");

	clen <- read.table(filename, header=T);
	mm <- simpleStd(clen[[1]], clen[[2]]);
	plen <- read.table(pfilename, header=T);
	plot(clen[[1]], clen[[2]], type="l", xlab="Length", ylab="Count");
	points(plen[[1]], plen[[2]], type="p");
	text(70, max(clen[[2]])*0.9, paste("mean", mm[1], sep="="))
	text(70, max(clen[[2]])*0.8, paste("std", mm[2], sep="="))
	text(70, max(clen[[2]])*0.7, paste("n", mm[3], sep="="))
	lengthsave <- c(lengthsave, clen)
}

plot(lengthsave[[1]], lengthsave[[2]]/sum(lengthsave[[2]]), type='l', xlim=c(43,53), xlab="Complete Length", ylab="Frequency");
points(lengthsave[[3]], lengthsave[[4]]/sum(lengthsave[[4]]), type='p', pch=2, col="red");
points(lengthsave[[5]], lengthsave[[6]]/sum(lengthsave[[6]]), type='p', pch=3, col="green");

dirs=c("shRNA_2_2_8ul_044", "shRNA_T5_protocol_run_047", "user_xh4_reanalyze_045");
labels <- c("228_044", "T5_047", "xh4_045")
-- plot complete
par(mfrow=c(2,2))
xlimit <- c(43,53)
i <- 1
for (d in dirs) {
	filename <- paste(d, "insertlencomplete.tab", sep="/");
	clen <- read.table(filename, header=T);
	mm <- dumStd(clen);
	topc <- max(clen[[2]])
	plot(clen[[1]], clen[[2]], type="l", xlab="Length", ylab="Count", xlim=xlimit);
	text(mean(xlimit)+2, topc*0.9, paste("mean", mm[1], sep="="))
	text(mean(xlimit)+2, topc*0.8, paste("std", mm[2], sep="="))
	text(mean(xlimit)+2, topc*0.7, paste("n", mm[3], sep="="))
	text(xlimit[1] + 2, topc*0.9, labels[i])
	i <- i+1
}


--point(lengthsave[[3]], lengthsave[[4]], type='p', lty=2);

-- for new runs we have different directory structure.
/*
treated:
        a118
		  a228
		  orig

untreated:
        a118
		  a228
		  orig
*/

dumMean <- function(data) {
	sum(data[[1]]*data[[2]])/sum(data[[2]])
}

dumStd <- function(data) {
	avg <- dumMean(data);
	total <- sum(data[[2]]);
	std <- sqrt(sum((data[[2]]/total)*(data[[1]] - avg)**2))
	c(round(avg, 3), round(std, 3), total)
}

dirs=c("untreated", "treated");
subdirs=c("orig", "a118", "a228");

-- complete sequences
xlimit <- c(45, 51)
par(mfrow=c(2,3))
for (d in dirs) {
	for (sd in subdirs) {
		filepath=paste(d, sd, "insertlencomplete.tab", sep="/");
		clen <- read.table(filepath, header=T);
		plot(clen[[1]], clen[[2]], type="l", xlab="Length", ylab="Count", xlim=xlimit);
		text(mean(xlimit) + 2, max(clen[[2]])*0.9, paste(d, sd, sep="-"))
		mm <- dumStd(clen)
		text(mean(xlimit) + 2, max(clen[[2]])*0.8, paste("mean", mm[1], sep="="))
		text(mean(xlimit) + 2, max(clen[[2]])*0.75, paste("std", mm[2], sep="="))
		text(mean(xlimit) + 2, max(clen[[2]])*0.7, paste("n", mm[3], sep="="))
	}
}

-- partial sequences
xlimit <- c(0, 100)
par(mfrow=c(2,3))
for (d in dirs) {
	for (sd in subdirs) {
		filepath=paste(d, sd, "insertlenpartial.tab", sep="/");
		plen <- read.table(filepath, header=T);
		mm <- dumStd(plen)
		plot(plen[[1]], plen[[2]], type="l", xlab="Length", ylab="Count", xlim=xlimit);
		text(xlimit[2]*0.65, max(plen[[2]])*0.8, paste(d, sd, sep="-"), adj=1, pos=3)
		text(mean(xlimit) + 2, max(plen[[2]])*0.8, paste("mean", mm[1], sep="="))
		text(mean(xlimit) + 2, max(plen[[2]])*0.75, paste("std", mm[2], sep="="))
		text(mean(xlimit) + 2, max(plen[[2]])*0.7, paste("n", mm[3], sep="="))
	}
}


