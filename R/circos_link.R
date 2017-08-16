#' A modified circos plot from OmicCircos package
#'
#' To draw circular plot with links shows contamination paired samples
#'
#' @param mapping, data frame
#' @param xc integer, circle center x coordinate
#' @param yc integer, circle center y coordinate
#' @param R interger, circle radius
#' @param W interger, circle width
#' @param cir, data frame
#' @param type character, plot type [chr|link]
#' @param print.chr.lab logical, TRUE: print chr label, FALSE: not print chr label
#' @param cex numeric, plot font size
#' @param lwd interger, line width
#' @param col vector, colors
#' @param side character, segment label position [in|out]
#' @param seg.lab.size numeric, segment label font size (default:1)
#' @param link.wd numeric, link width (default:1)
#' @param lower_arch numeric, (default:1)
#' @param line_type character, line type (default: "l")
#'
#' @usage circos_link <- function (mapping=mapping, xc=400, yc=400, R=400, W=W,cir="", type="n", print.chr.lab=TRUE, cex=1,lwd=1, col=rainbow(10, alpha=0.8)[7], side="",seg.lab.size=1, link.wd=1,lower_arch=1,line_type="l")
#' @return .tif figure at output directory
#'
#' @references
#' {TBA}
#'
#' @export
#'

## main ##
circos_link <- function (mapping=mapping, xc=400, yc=400, R=400, W=W,
    cir="", type="n", print.chr.lab=TRUE, cex=1,lwd=1, col=rainbow(10, alpha=0.8)[7], side="", scale=FALSE,
    seg.lab.size=1, link.wd=1,lower_arch=1,line_type="l"
  ){

  ############################ main #############################
  # initial
  r        <- R;

  #### cir ####
  # In cir, the segments (chromosomes) should be sorted.
  cir.m <- 0;

  if (is.matrix(cir)==T|is.data.frame(cir)==T|class(cir)[1] == "GRanges") {
      # based on the input file,
      if (class(cir)[1] == "GRanges"){
        chr.po    <- cbind(as.character(seqnames(cir)),
                    start(cir), names(cir), mcols(cir));
      } else {
        chr.po    <- cir;
      }
      chr.po[,1]  <- gsub("chr","",chr.po[,1]);
      chr.num     <- nrow(chr.po);
      cir.m       <- 1;
  } else if (nchar(cir) >= 1){
      # by known genome
      pofile.s <- paste("UCSC.", cir, ".RData", sep="");
      pofile   <- system.file("data", pofile.s, package="OmicCircos");
      if (file.exists(pofile)){
      } else {
        stop ("cir name?");
      }
      chr.po     <- local(get(load(pofile)));
      chr.po[,1] <- gsub("chr","",chr.po[,1]);
      chr.num    <- nrow(chr.po);
  } else {
      stop("cir name ??")
  }
  # end cir

  #### chr ###
  if (type == "chr"){
    chr.lw = W;
    if (cir.m == 0){
      # for chromosomes have cytoband data
      chrfile.s <- paste("UCSC.", cir, ".chr.RData", sep="");
      chrfile   <- system.file("data", chrfile.s, package="OmicCircos");
      if (file.exists(chrfile)){
      } else {
        stop ("chrfile name?");
      }
      dat.c     <- local(get(load(chrfile)));
      dat.c[,1] <- gsub("chr", "", dat.c[,1]);

      for (chr.i in c(1:chr.num)){
        chr.s  <- chr.po[chr.i,1];

        v1 <- as.numeric(chr.po[chr.i,2]);
        v2 <- as.numeric(chr.po[chr.i,3]);
        v3 <- as.numeric(chr.po[chr.i,6]);
        v4 <- as.numeric(chr.po[chr.i,7]);

        dat.v <- subset(dat.c, dat.c[,1]==chr.s);
        dat.v <- dat.v[order(as.numeric(dat.v[,2])),];
        for (i in 1:nrow(dat.v)){
          if (dat.v[i,5]=="gneg"){
            col <- colors()[351];
          } else if (dat.v[i,5]=="acen" | dat.v[i,5]=="gvar" | dat.v[i,5]=="stalk"){
            col <- colors()[26];
          } else {
            col.v <- gsub("gpos","",dat.v[i,5]);
            if (col.v >= 30){
              col <- colors()[300];
            } else {
              col <- colors()[351];
            }
          }
          w1 <- scale.v(dat.v[i,2], v1, v2, v3, v4);
          w2 <- scale.v(dat.v[i,3], v1, v2, v3, v4);
          draw.arc.s(xc, yc, r, w1, w2, col=col[i], lwd=chr.lw);
        }
        if (print.chr.lab){
          v1  <- as.numeric(chr.po[chr.i,2]);
          v2  <- as.numeric(chr.po[chr.i,3]);
          w.m <- (v1+v2)/2;
          r.j <- r/20;
          chr.t <- gsub("chr", "", chr.s);
          draw.text.rt(xc, yc, r+r.j, w.m, chr.t, cex=seg.lab.size); ## LP
        }
        if (scale){
          total.num <- as.numeric(chr.po[nrow(chr.po),5]);
          do.scale.cir(xc=xc, yc=yc, the.r=r+10, total.num=total.num,
                       col="blue", lwd=0.001,
                       V1=v1, V2=v2, V3=v3, V4=v4);
        }
      }
    } else {
        if (!is.null(col)){
          if(length(col)==1){
          col <- rep(col, nrow(chr.po))[c(1:nrow(chr.po))];
          }
        } else {
          col <- rep("blue", nrow(chr.po));
        }

        for (chr.i in c(1:chr.num)){
          w1 <- as.numeric(chr.po[chr.i,2]);
          w2 <- as.numeric(chr.po[chr.i,3]);
          draw.arc.s(xc, yc, r, w1, w2, col=as.character(col[chr.i]), lwd=chr.lw); ##LP

          if (print.chr.lab){
            w.m <- (w1+w2)/2;
            r.j <- W/2 + r/10;
            chr.t <- gsub("chr", "", chr.po[chr.i,1]);
            draw.text.rt(xc, yc, r+r.j, w.m, chr.t, cex=seg.lab.size); ## LP
          }

          if (scale){
            v1 <- as.numeric(chr.po[chr.i,2]);
            v2 <- as.numeric(chr.po[chr.i,3]);
            v3 <- as.numeric(chr.po[chr.i,6]);
            v4 <- as.numeric(chr.po[chr.i,7]);

            total.num <- as.numeric(chr.po[nrow(chr.po),5]);
            do.scale.cir(xc=xc, yc=yc, the.r=r+10, total.num=total.num,
                         col="blue", lwd=0.001,
                         V1=v1, V2=v2, V3=v3, V4=v4);
          }
       }
    }
  }
  ### end chr ###

  ### type  ################


   ### link
  if (type == "link"){
    chr.po[,4] <- gsub("chr", "", chr.po[,4]);
    dat.in <- mapping;
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    dat.in[,4] <- gsub("chr", "", dat.in[,4]);
    dat    <- dat.in;

    if (!is.null(col)){
      col <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];
    }
    if (!is.null(lwd)){
      ##lwd <- rep(lwd, nrow(dat.in))[c(1:nrow(dat.in))];
      lwd <-round(as.numeric(link.wd)+1,digits=0) ##LP

    } else {
      lwd <- rep(1, nrow(dat.in));
    }

    r <- R;
    #checkn= nrow(dat)
    for (i in 1:nrow(dat)){
      chr1.s   <- dat[i,1];
      chr2.s   <- dat[i,4];
      po1      <- dat[i,2];
      po2      <- dat[i,5];

      chr1     <- which(chr.po[,1]==chr1.s);
      chr2     <- which(chr.po[,1]==chr2.s);
      #print(c(i,checkn,chr1,chr2))
      v1 <- as.numeric(chr.po[chr1,2]);
      v2 <- as.numeric(chr.po[chr1,3]);
      v3 <- as.numeric(chr.po[chr1,6]);
      v4 <- as.numeric(chr.po[chr1,7]);

      w1 <- scale.v(as.numeric(po1), v1, v2, v3, v4);

      v1 <- as.numeric(chr.po[chr2,2]);
      v2 <- as.numeric(chr.po[chr2,3]);
      v3 <- as.numeric(chr.po[chr2,6]);
      v4 <- as.numeric(chr.po[chr2,7]);

      w2 <- scale.v(as.numeric(po2), v1, v2, v3, v4);

      if (chr1 == chr2){
        draw.link(xc, yc, r, w1, w2, col=as.character(col[i]), lwd=lwd[i],lower_arch = lower_arch,line_type = line_type); #LP
      } else {
        r.j <- W/2; ##LP
        draw.link(xc, yc, r-r.j, w1, w2, col=as.character(col[i]), lwd=lwd[i],lower_arch=lower_arch ,line_type=line_type[i]); ##LP

      }
    }
  }
  ### end link
  if (type == "highlight.link"){
     xc   <- mapping[1];
     yc   <- mapping[2];
     r    <- mapping[3];
     w1.1 <- mapping[4];
     w1.2 <- mapping[5];
     w2.1 <- mapping[6];
     w2.2 <- mapping[7];
     draw.link.pg(xc, yc, r, w1.1, w1.2, w2.1, w2.2, col=col, lwd=lwd);
  }

### end type  ################
}

##############################
### functions
##############################
bezierCurve <- function(x, y, n=10)	{
  outx <- NULL
  outy <- NULL
  i <- 1

  for (t in seq(0, 1, length.out=n))		{
	b <- bez(x, y, t)
    outx[i] <- b$x
    outy[i] <- b$y
    #print(c(i,t))
	  i <- i+1
  }
  #print(cbind(outx,outy))
  return (list(x=outx, y=outy))
}

###
bez <- function(x, y, t)	{
  outx <- 0
  outy <- 0
  n <- length(x)-1
  #print(c("t:",t,"n:",n,"x:",x,"y:",y))
  for (i in 0:n)		{
    outx <- outx +  choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]
	  outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]
	 }
  return (list(x=outx, y=outy))
}

###########################################
# one value : from a to b
scale.v <- function(v, a, b, min.v, max.v) {
  v <- v-min.v;
  v <- v/(max.v-min.v);
  v <- v*(b-a);
  v+a
}

## permutation index
perm_list <- function (n, r, v = 1:n){
    if (r == 1)
       X <- matrix(v, n, 1)
    else if (n == 1)
       X <- matrix(v, 1, r)
    else {
       X <- NULL
       for (i in 1:n){
            X <- rbind(X, cbind(v[i], perm_list(n-1 , r-1 , v[-i])))
       }
    }
    return(X);
}

### draw.link
draw.link <- function(xc, yc, r, w1, w2, col=col, lwd=lwd,lower_arch,line_type='l') {
    # for translocation
    w3  <- (w1+w2)/2;

    w1  <- w1/360*2*pi;
    w2  <- w2/360*2*pi;
    w3  <- w3/360*2*pi;

    x0  <- xc+r*cos(w1);
    y0  <- yc-r*sin(w1);
    x1  <- xc+r*cos(w2);
    y1  <- yc-r*sin(w2);

    x2  <- xc+(r/lower_arch)*cos(w3);##LP
    y2  <- yc-(r/lower_arch)*sin(w3);##LP

    if(lower_arch!=1){
      x <- c(x0,x2,x2,x1);##LP
      y <- c(y0,y2,y2,y1);##LP
    }else{
      x <- c(x0,xc,xc,x1);
      y <- c(y0,yc,yc,y1);
    }

      #print(x)
      #print(y)

    points(bezierCurve(x,y,20), col=col, lwd=lwd, lend="butt",type=as.character(line_type))
}

###
draw.link.pg <- function(xc, yc, r, w1.1, w1.2, w2.1, w2.2, col=col, lwd=lwd) {
    #####################################################
    w1 <- w1.1;
    w2 <- w2.2;
    w3  <- (w1+w2)/2;
    w1  <- w1/360*2*pi;
    w2  <- w2/360*2*pi;
    w3  <- w3/360*2*pi;
    x0  <- xc+r*cos(w1);
    y0  <- yc-r*sin(w1);
    x1  <- xc+r*cos(w2);
    y1  <- yc-r*sin(w2);
    x <- c(x0,xc,xc,x1);
    y <- c(y0,yc,yc,y1);
    bc1 <- bezierCurve(x,y,60);

    ang.d <- abs(w1.1-w1.2);
    pix.n <- ang.d * 10;
    if (pix.n < 10){
      pix.n <- 10;
    }

    ang.seq <- rev(seq(w1.1,w1.2,length.out=pix.n));
    ang.seq <- ang.seq/360*2*pi;

    fan.1.x <- xc + cos(ang.seq) * r;
    fan.1.y <- yc - sin(ang.seq) * r;

    ######################################################
    w1 <- w1.2;
    w2 <- w2.1;
    w3  <- (w1+w2)/2;
    w1  <- w1/360*2*pi;
    w2  <- w2/360*2*pi;
    w3  <- w3/360*2*pi;
    x0  <- xc+r*cos(w1);
    y0  <- yc-r*sin(w1);
    x1  <- xc+r*cos(w2);
    y1  <- yc-r*sin(w2);
    x <- c(x0,xc,xc,x1);
    y <- c(y0,yc,yc,y1);
    bc2 <- bezierCurve(x,y,60);

    ang.d <- abs(w2.1-w2.2);
    pix.n <- ang.d * 10;
    if (pix.n < 10){
      pix.n <- 10;
    }

    ang.seq <- rev(seq(w2.1,w2.2,length.out=pix.n));
    ang.seq <- ang.seq/360*2*pi;

    fan.2.x <- xc + cos(ang.seq) * r;
    fan.2.y <- yc - sin(ang.seq) * r;

    polygon(c(bc1$x, fan.2.x, rev(bc2$x), rev(fan.1.x)),
            c(bc1$y, fan.2.y, rev(bc2$y), rev(fan.1.y)),
            fillOddEven=TRUE, border=col, col=col, lwd=lwd);
}

###
draw.text.rt <- function(xc, yc, r, w, n, col="black", cex=1, side="out"){

  w     <- w%%360;
  the.o <- w;

  if (w <= 90){
    #w=0.0219780219780059+0.978021978021979*w;
  } else if (w > 90 & w <= 180){
    #w=3.91304347826081+0.978260869565218*w;
  } else if (w > 180 & w <= 270){
    #w=3.91304347826091+0.978260869565217*w;
  } else if (w > 270 & w <= 360){
    #w=7.82608695652158+0.978260869565218*w;
  }

  the.w <- 360-w;
  w     <- w/360*2*pi;
  x     <- xc+r*cos(w);
  y     <- yc-r*sin(w);

  num2  <- 26;

  if (side=="out"){
    if (the.w <= 90 ){
      the.pos <- 4;
    } else if (the.w > 90 & the.w <= 180) {
      the.w <- the.w + 180;
      the.pos <- 2;
    } else if (the.w > 180 & the.w <= 270){
      the.w <- the.w%%180;
      the.pos <- 2;
    } else if (the.w > 270 & the.w <= 360){
      the.pos <- 4;
    }

    if (the.pos==2){
      x <- x+num2;
    }
    if (the.pos==4){
      x <- x-num2;
    }
  }

  if (side=="in"){
    if (the.w <= 90 ){
      the.pos <- 4;
    } else if (the.w > 90 & the.w <= 180) {
      the.w <- the.w + 180;
      the.pos <- 2;
    } else if (the.w > 180 & the.w <= 270){
      the.w <- the.w%%180;
      the.pos <- 2;
    } else if (the.w > 270 & the.w <= 360){
      the.pos <- 4;
    }

    if (the.pos==2){
      x <- x+num2;
    }
    if (the.pos==4){
      x <- x-num2;
    }
  }

  text(x, y, adj=0, offset=1, labels=n, srt=the.w,
       pos=the.pos, col=col, cex=cex);
}

###strokeLine2
draw.line <- function (xc, yc, w, l1, l2, col=col, lwd=lwd, lend=1) {
    w  <- (w/360)*2*pi;
    x1 <- xc+l1*cos(w);
    y1 <- yc-l1*sin(w);
    x2 <- xc+l2*cos(w);
    y2 <- yc-l2*sin(w);
    segments(x1, y1, x2, y2, col=col, lwd=lwd, lend=lend);
}

###strokeLine3
draw.line2 <- function (xc, yc, w, r, l, col=col, lwd=lwd){
    line_w   <- l;
    theangle <- w;
    l1       <- r;
    theangle <- (theangle/360)*2*pi;
    x0       <- xc+l1*cos(theangle);
    y0       <- yc+l1*sin(theangle);
    w1       <- 45/360*2*pi;
    x1 = xc + sin(w1) * (x0);
    y1 = yc + cos(w1) * (y0);
    x2 = xc - sin(w1) * (x0);
    y2 = yc - cos(w1) * (y0);
    segments(x1, y1, x2, y2, col=col, lwd=lwd, lend="butt");
}

###strokeLine by two angles
draw.line3 <- function (xc, yc, w1, w2, r1, r2, col=col, lwd=lwd){
    theangle1 <- w1;
    theangle2 <- w2;
    l1        <- r1;
    l2        <- r2;

    theangle1 <- (theangle1/360)*2*pi;
    x1        <- xc+l1*cos(theangle1);
    y1        <- yc-l1*sin(theangle1);

    theangle2 <- (theangle2/360)*2*pi;
    x2        <- xc+l2*cos(theangle2);
    y2        <- yc-l2*sin(theangle2);

    segments(x1, y1, x2, y2, col=col, lwd=lwd, lend="butt");
}

draw.arc.s <- function (xc, yc, r, w1, w2, col="lightblue", lwd=1, lend=1){
  # s = simple
  # r = radius

  ang.d <- abs(w1-w2);
  pix.n <- ang.d * 5;
  if (pix.n < 2){
    pix.n <- 2;
  }

  ang.seq <- rev(seq(w1,w2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;

  fan.i.x <- xc + cos(ang.seq) * r;
  fan.i.y <- yc - sin(ang.seq) * r;
  ## lend=0(round); lend=1(butt); lend=2(square)
  lines(fan.i.x, fan.i.y, col=col, lwd=lwd, type="l", lend=lend);
  #points(fan.i.x, fan.i.y, col=col, lwd=lwd, type="l", lend=lend);
  #print(paste("plot:",col))
}
#########################################
## segment to angle and position
## segAnglePo
#########################################

# get angle if given seg and position
# seg=chromosome in genome, shuch as: 1:22, 23="X", 24="Y";

setGeneric("segAnglePo", function(seg.dat=seg.dat, seg=seg, angle.start=angle.start,
                                  angle.end=angle.end){standardGeneric("segAnglePo")})

setMethod("segAnglePo", "data.frame", function(seg.dat=seg.dat, seg=seg,
                                               angle.start=angle.start, angle.end=angle.end){
	  ##
	  seg.dat  <- seg.dat;
	 .segAnglePo(seg.dat=seg.dat, seg=seg, angle.start=angle.start, angle.end=angle.end);
})

setMethod("segAnglePo", "GRanges", function(seg.dat=seg.dat){
	  ##
	  seg.dat  <- cbind(as.character(seqnames(seg.dat)), start(seg.dat), end(seg.dat),
                      names(seg.dat), as.character(strand(seg.dat)));
	  colnames(seg.dat) <- c("seg.name","seg.Start","seg.End","name","gieStain");
	  .segAnglePo(seg.dat=seg.dat, seg=seg, angle.start=angle.start, angle.end=angle.end);
})

## seg should be ordered by user
.segAnglePo <- function (seg.dat=seg.dat, seg=seg, angle.start=angle.start,
                         angle.end=angle.end){

  if (missing(angle.start)){
    angle.start <- 0;
  }
  if (missing(angle.end)){
    angle.end <- 360;
  }
## check data.frame?
  colnames(seg.dat) <- c("seg.name","seg.Start","seg.End","name","gieStain");

## get length of the segomosomes
  seg.l   <- c();
  seg.min <- c();
  seg.sum <- c();
  seg.s   <- 0;
  seg.num   <- length(seg);
  seg.names <- seg;

########################################################
########################################################
  for (i in 1:seg.num){
    seg.n <- seg.names[i];
    dat.m <- subset(seg.dat, seg.dat[,1]==seg.n);
    seg.full.l <- max(as.numeric(dat.m[,"seg.End"]));
    seg.full.m <- min(as.numeric(dat.m[,"seg.Start"]));
    seg.l      <- cbind(seg.l, seg.full.l);
    seg.min    <- cbind(seg.min, seg.full.m);
    seg.s      <- seg.s + seg.full.l;
    seg.sum    <- cbind(seg.sum, seg.s);
 }

 ## initial parameters
 gap.angle.size <- 2;
 seg.angle.from <- angle.start + 270;
 seg.full.l  <- sum(as.numeric(seg.l));
 angle.range <- angle.end - angle.start;
 cir.angle.r <- (angle.range - seg.num * gap.angle.size)/seg.full.l;
 out.s     <- c();
 l.old     <- 0;
 gap.angle <- 0;

 for (i in 1:seg.num){
   seg.n <- seg.names[i];
   dat.m <- subset(seg.dat, seg.dat[,1]==seg.n);
   len   <- seg.sum[i];
   w1    <- cir.angle.r*l.old + gap.angle;
   w2    <- cir.angle.r*len   + gap.angle;
   out.s     <- rbind(out.s, c(seg.n, w1+seg.angle.from, w2+seg.angle.from, l.old, len, seg.min[i], seg.l[i]));
   gap.angle <- gap.angle + gap.angle.size;
   l.old     <- len;
 }
  return(out.s);
}

sim.circos <- function (seg=10, po=c(20,50), ind=10, link=10,
        link.pg=10){
 seg.num    <- seg;
 seg.po     <- po;
 ind.num    <- ind;
 link.num   <- link;
 link.w.num <- link.pg;

 ################################################
 ## segment, pointer, individual
 ################################################
 seg.out <- c();
 seg.v   <- c();

 for (i in 1:seg.num){
  po.num <- sample(seg.po, 1);
  i.s    <- paste("chr", i, sep="");
  for (j in 1:po.num){
    seg.out <- rbind(seg.out, c(i.s, j-1, j, "NA", "NA"));
    c.v <- c();
    for (k in 1:ind.num){
      if (k > ind.num/2){
        the.v   <- round(rnorm(1) + i + i/ind.num , 3);
      } else {
        the.v   <- round(rnorm(1) + i,3);
      }
      c.v     <- c(c.v, the.v);
    }
    seg.v   <- rbind(seg.v,   c(i.s, j, c.v));
  }
 }

 colnames(seg.out) <- c("seg.name", "seg.Start", "seg.End", "the.v", "NO")
 seg.out <- as.data.frame(seg.out);
 names   <- paste("name", 1:ind.num, sep="");
 colnames(seg.v) <- c("seg.name", "seg.po", names)
 seg.v   <- as.data.frame(seg.v);

 ################################################
 ## link data
 ################################################
 link.out <- c();
 for (i in 1:link.num){
  name <- paste("n", i, sep="");
  c1 <- sample(c(1:nrow(seg.out)), 1);
  c2 <- sample(c(1:nrow(seg.out)), 1);
  link.out <- rbind(link.out, c(seg.out[c1,1], seg.out[c1,2], name,
                              seg.out[c2,1], seg.out[c2,2], name, name));
 }
 colnames(link.out) <- c("seg1", "po1", "name1", "seg2", "po2", "name2", "name3")
 link.out <- as.data.frame(link.out);

 ################################################
 ## link data with wide
 ################################################
 link.w.out <- c();
 for (i in 1:link.w.num){
  n.i <- sample(1:seg.num, 2);
  i1    <- paste("chr", n.i[1], sep="");
  i2    <- paste("chr", n.i[2], sep="");
  n1    <- subset(seg.out, seg.out[,1]==i1);
  n2    <- subset(seg.out, seg.out[,1]==i2);
  re1   <- sample(n1[,2], 2);
  re2   <- sample(n2[,2], 2);
  link.w.out <- rbind(link.w.out, c(i1, re1, i2, re2));
 }
 colnames(link.w.out) <- c("seg1", "start1", "end1", "seg2", "start2", "end2");
 link.w.out <- as.data.frame(link.w.out)

 sim.out <- list(seg.frame=seg.out, seg.mapping=seg.v,
    seg.link = link.out, seg.link.pg = link.w.out);

 return(sim.out);
}


## color bar
## xl=xleft, yb=ybottom, xr=xright, yt=ytop
color.bar <- function(xl, yb, xr, yt, v.min, v.max) {
    nticks <- 11;
    ticks  <- seq(v.min, v.max, len=nticks);

    lut   <- colorRampPalette(c("blue", "white", "red"))(100);
    scale <- (length(lut)-1)/(v.max-v.min);
    ys    <- (yt-yb)/(length(lut)-1);

    yb.old <- yb;
    for (i in 1:(length(lut)-1)) {
      rect(xl, yb.old, xr, yb.old+ys, col=lut[i], border=NA)
      yb.old <- yb.old + ys;
    }
    v.med <- (v.max+v.min)/2;

    ym  <- (yt+yb)/2;

    v.max <- round(v.max,2);
    v.min <- round(v.min,2);
    v.med <- round(v.med,2);
    text(xl, yt+30, "heatmap", cex=1);
    text(xl-30, yt, v.max, cex=0.6);
    text(xl-30, yb, v.min, cex=0.6);
    text(xl-30, ym, v.med, cex=0.6);
    segments(xl-5, yb, xl-5,  yt, col="black", lwd=0.8);
    segments(xl-5, yb, xl-15, yb, col="black", lwd=0.8);
    segments(xl-5, yt, xl-15, yt, col="black", lwd=0.8);
    segments(xl-5, ym, xl-15, ym, col="black", lwd=0.8);

}

### heatmap.cluster
heatmap.cluster <- function(x1, y1, x2, y2, dat=dat){
    dat.d  <- dist(t(dat));
    dat.h  <- hclust(dat.d);

    ####################
    c.x1   <- x1;
    c.x2   <- x2;
    c.y1   <- y1;
    c.y2   <- y2;
    max.h  <- max(dat.h$height);
    max.n  <- length(dat.h$labels);

    lab   = (c.x2-c.x1)/max.n;
    ratio = (c.y2-c.y1)/max.h;

    ####################
    sub2h  <- c();
    sub2po <- c();
    y0 <- c.y1 - 2;
    for (i in 1:nrow(dat.h$merge)){
      le <- dat.h$merge[i,1];
      ri <- dat.h$merge[i,2];
      if (le < 0 && ri < 0){
        po1 <- which(dat.h$order==abs(le));
        po2 <- which(dat.h$order==abs(ri));
        y1  <- y0;
        y2  <- y0;
      } else if (le < 0){
        po1 <- which(dat.h$order==abs(le));
	   po2 <- sub2po[ri,2];
	   y1  <- y0;
	   y2  <- c.y1 + ratio * sub2h[ri,2];
      } else if (ri < 0){
	   po1 <- sub2po[le, 2];
	   po2 <- which(dat.h$order==abs(ri));
	   y1  <- c.y1 + ratio * sub2h[le,2];
	   y2  <- y0;
      } else {
	   po1 = sub2po[le, 2];
	   po2 = sub2po[ri, 2];
	   y1  = c.y1 + ratio * sub2h[le,2];
	   y2  = c.y1 + ratio * sub2h[ri,2];
      }
      sub2po <- rbind(sub2po, c(i, (po1+po2)/2));
      sub2h  <- rbind(sub2h,  c(i, dat.h$height[i]));
      x1 = c.x1 + lab*po1;
      x2 = c.x1 + lab*po2;
      y3 = c.y1 + ratio*dat.h$height[i];
      segments(x1,y3,x2,y3, lwd=0.1);
      # le
      segments(x1,y3,x1,y1, lwd=0.1);
      # ri
      segments(x2,y3,x2,y2, lwd=0.1);
    }

    y0 <- c.y1 - 8;
    for (i in 1:length(dat.h$order)){
       n  <- dat.h$label[dat.h$order[i]];
       x0 <- c.x1 + lab * i - 24;
       text(x0, y0, n, srt=270, col="blue", cex=0.2, pos=4, adj=0, offset=1);
    }
    ####################
}
### end heatmap.cluster

### start do.scale
do.scale <- function(xc=xc, yc=yc, dat.min=dat.min, dat.max=dat.max,
   R=R, W=W, s.n=1, col="blue"){
  dat.m   <- round((dat.min+dat.max)/2, s.n);
  dat.min <- round(dat.min, s.n);
  dat.max <- round(dat.max, s.n);
  y1      <- yc + R ;
  y2      <- yc + R + W/2;
  y3      <- yc + R + W;
  x1      <- xc - W/20;
  x2      <- x1 - (W/20)*1.2;
  x3      <- x1 - (W/20)*3;
  segments(x1, y1, x1, y3, lwd=0.01, col=col);
  segments(x1, y1, x2, y1, lwd=0.01, col=col);
  segments(x1, y2, x2, y2, lwd=0.01, col=col);
  segments(x1, y3, x2, y3, lwd=0.01, col=col);
  text(x3, y1, dat.min, cex=0.2, col=col);
  text(x3, y2, dat.m,   cex=0.2, col=col);
  text(x3, y3, dat.max, cex=0.2, col=col);
}
### end do.scale

### start do.scale.cir
do.scale.cir <- function (xc=xc, yc=yc, the.r=the.r, total.num=total.num,
                col="blue", lwd=0.001,
                V1=V1, V2=V2, V3=V3, V4=V4
             ){

  ### scale initial start ##################
  ## scale label number = 150 in 360 degree circumference
  ## the.r  = R
  ## sum.po = total point number
  ## should return: 1) scale number 2) scale unit
  sum.po  <- as.numeric(total.num);         # total number
  scale.w <- sum.po/150;                         # one scale size
  scale.l <- nchar(as.integer(scale.w))-1;       # length of the number
  scale.i <- as.integer(scale.w/(10^scale.l));   # the first digital
  scale.d <- scale.i * 10^scale.l;               # the first digital, then 0
  scale.m <- 1 * 10^scale.l;                     # the unit of the scale
  ### scale initial end   ##################

  ### scale start ######################
  draw.arc.s(xc, yc, w1=V1, w2=V2, r=the.r, col=col, lwd=lwd);

  start.p  <- 0;
  w        <- scale.v(as.numeric(start.p), V1, V2, V3, V4);
  scale.n  <- as.integer(as.numeric(V4)/scale.d+1);

  for (j in 1:scale.n){
    po.s <- as.integer(start.p/scale.m);

    draw.line(xc,    yc,   w, the.r, the.r+4, col=col, lwd=lwd);
    draw.text.rt(xc, yc,  the.r+6, w, po.s, col=col, cex=0.2);

    start.p <- start.p + scale.d/2;
    w        <- scale.v(as.numeric(start.p), V1, V2, V3, V4);

    if (w <= V2){
      draw.line(xc, yc, w, the.r, the.r+2, col=col, lwd=lwd);
    }

    start.p <- start.p + scale.d/2;
    w        <- scale.v(as.numeric(start.p), V1, V2, V3, V4);
  }
  ### scale end ######################
}
### end do.scale.cir

### zoom ###
zoom.in <- function(cir.in=cir.in, zoom=zoom){
   out.cir  <- c();
   out.dat  <- c();

   ## for all data sets
   chr1     <- which(cir.in[,1]==zoom[1]);
   chr2     <- which(cir.in[,1]==zoom[2]);

   v1 <- as.numeric(cir.in[chr1,2]);
   v2 <- as.numeric(cir.in[chr1,3]);
   v3 <- as.numeric(cir.in[chr1,6]);
   v4 <- as.numeric(cir.in[chr1,7]);
   w3 <- scale.v(zoom[3], v1, v2, v3, v4);

   v1 <- as.numeric(cir.in[chr2,2]);
   v2 <- as.numeric(cir.in[chr2,3]);
   v3 <- as.numeric(cir.in[chr2,6]);
   v4 <- as.numeric(cir.in[chr2,7]);
   w4 <- scale.v(zoom[4], v1, v2, v3, v4);

   for (i in chr1:chr2){
     w1 <- as.numeric(zoom[5])+270;
     w2 <- as.numeric(zoom[6])+270;

     v1 <- scale.v(as.numeric(cir.in[i,2]), w1, w2, w3, w4);
     v2 <- scale.v(as.numeric(cir.in[i,3]), w1, w2, w3, w4);
     out.cir <- rbind(out.cir, c(as.character(cir.in[i,1]), v1, v2, cir.in[i,4:ncol(cir.in)]));
   }
   return(out.cir);
}
### end zoom ###

###################################################################
## from gtools package
mixedorder <- function (x)
{
    delim = "\\$\\@\\$"
    numeric <- function(x) {
        optwarn = options("warn")
        on.exit(options(optwarn))
        options(warn = -1)
        as.numeric(x)
    }
    nonnumeric <- function(x) {
        optwarn = options("warn")
        on.exit(options(optwarn))
        options(warn = -1)
        ifelse(is.na(as.numeric(x)), toupper(x), NA)
    }
    x <- as.character(x)
    which.nas <- which(is.na(x))
    which.blanks <- which(x == "")
    if (length(which.blanks) > 0)
        x[which.blanks] <- -Inf
    if (length(which.nas) > 0)
        x[which.nas] <- Inf
    delimited <- gsub("([+-]{0,1}[0-9]+([eE][\\+\\-]{0,1}[0-9]+){0,1})",
        paste(delim, "\\1", delim, sep = ""), x)
    step1 <- strsplit(delimited, delim)
    step1 <- lapply(step1, function(x) x[x > ""])
    step1.numeric <- lapply(step1, numeric)
    step1.character <- lapply(step1, nonnumeric)
    maxelem <- max(sapply(step1, length))
    step1.numeric.t <- lapply(1:maxelem, function(i) sapply(step1.numeric,
        function(x) x[i]))
    step1.character.t <- lapply(1:maxelem, function(i) sapply(step1.character,
        function(x) x[i]))
    rank.numeric <- sapply(step1.numeric.t, rank)
    rank.character <- sapply(step1.character.t, function(x) as.numeric(factor(x)))
    rank.numeric[!is.na(rank.character)] <- 0
    rank.character <- t(t(rank.character) + apply(matrix(rank.numeric),
        2, max, na.rm = TRUE))
    rank.overall <- ifelse(is.na(rank.character), rank.numeric,
        rank.character)
    order.frame <- as.data.frame(rank.overall)
    if (length(which.nas) > 0)
        order.frame[which.nas, ] <- Inf
    retval <- do.call("order", order.frame)
    return(retval)
}
