#-*- perl -*-

# Copyright (C) 2001--2007 R Development Core Team
# Copyright (C) 2003-4, 2006 The R Foundation
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# A copy of the GNU General Public License is available at
# http://www.r-project.org/Licenses/
#
# Send any bug reports to r-bugs@r-project.org

## Usage: perl massage-Examples.pl pkgname files

## Given a list of files of the form .../.../<name>.R, produce one large
## file, i.e., write to stdout, concatenating the files together with
## 1) Putting a HEADER in front
## 2) Wrapping every file in order to be more order independent
## 3) appending a FOOTER ...

use File::Basename;

my $PKG = shift @ARGV;
my @Rfiles;
if(-d $ARGV[0]) {
    my $dir = $ARGV[0];
    opendir(DIR, $dir) or die "cannot opendir $dir: $!";
    my @files = sort grep { /\.R$/ } readdir(DIR);
    closedir(DIR);
    foreach my $file (@files) {
	push(@Rfiles, "$dir/$file");
    }
} else {
    @Rfiles = @ARGV;
}

### * Header
print <<_EOF_;
### * <HEADER>
###
attach(NULL, name = "CheckExEnv")
assign("nameEx", 
       local({
	   s <- "__{must remake R-ex/*.R}__"
           function(new) {
               if(!missing(new)) s <<- new else s
           }
       }),
       pos = "CheckExEnv")
## Add some hooks to label plot pages for base and grid graphics
assign("base_plot_hook",
       function() {
           pp <- par(c("mfg","mfcol","oma","mar"))
           if(all(pp\$mfg[1:2] == c(1, pp\$mfcol[2]))) {
               outer <- (oma4 <- pp\$oma[4]) > 0; mar4 <- pp\$mar[4]
               mtext(sprintf("help(\\"%s\\")", nameEx()), side = 4,
                     line = if(outer)max(1, oma4 - 1) else min(1, mar4 - 1),
               outer = outer, adj = 1, cex = .8, col = "orchid", las=3)
           }
       },
       pos = "CheckExEnv")
assign("grid_plot_hook",
       function() {
           grid::pushViewport(grid::viewport(width=grid::unit(1, "npc") - 
                              grid::unit(1, "lines"), x=0, just="left"))
           grid::grid.text(sprintf("help(\\"%s\\")", nameEx()),
                           x=grid::unit(1, "npc") + grid::unit(0.5, "lines"),
                           y=grid::unit(0.8, "npc"), rot=90,
                           gp=grid::gpar(col="orchid"))
       },
       pos = "CheckExEnv")
setHook("plot.new",     get("base_plot_hook", pos = "CheckExEnv"))
setHook("persp",        get("base_plot_hook", pos = "CheckExEnv"))
setHook("grid.newpage", get("grid_plot_hook", pos = "CheckExEnv"))
assign("cleanEx",
       function(env = .GlobalEnv) {
	   rm(list = ls(envir = env, all.names = TRUE), envir = env)
           RNGkind("default", "default")
	   set.seed(1)
   	   options(warn = 1)
	   .CheckExEnv <- as.environment("CheckExEnv")
	   delayedAssign("T", stop("T used instead of TRUE"),
		  assign.env = .CheckExEnv)
	   delayedAssign("F", stop("F used instead of FALSE"),
		  assign.env = .CheckExEnv)
	   sch <- search()
	   newitems <- sch[! sch %in% .oldSearch]
	   for(item in rev(newitems))
               eval(substitute(detach(item), list(item=item)))
	   missitems <- .oldSearch[! .oldSearch %in% sch]
	   if(length(missitems))
	       warning("items ", paste(missitems, collapse=", "),
		       " have been removed from the search path")
       },
       pos = "CheckExEnv")
assign("ptime", proc.time(), pos = "CheckExEnv")
## at least one package changes these via ps.options(), so do this
## before loading the package.
## Use postscript as incomplete files may be viewable, unlike PDF.
## Choose a size that is close to on-screen devices, fix paper
ps.options(width = 7, height = 7, paper = "a4", reset = TRUE)
grDevices::postscript("$PKG-Ex.ps")
		      
assign("par.postscript", graphics::par(no.readonly = TRUE), pos = "CheckExEnv")
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
options(warn = 1)    
_EOF_

if ($PKG eq "tcltk") {
    print "require('tcltk') || q()\n\n";
} elsif ($PKG ne "base") {
    print "library('$PKG')\n\n";
}
print "assign(\".oldSearch\", search(), pos = 'CheckExEnv')\n";
print "assign(\".oldNS\", loadedNamespaces(), pos = 'CheckExEnv')\n";

### * Loop over all R files, and edit a few of them ...
foreach my $file (@Rfiles) {
    my $have_examples = 0;
    my $have_par = 0;
    my $have_contrasts = 0;
    my $nm;

    $nm = basename $file, (".R");
    $nm =~ s/[^- .a-zA-Z0-9]/./g;

    if ($PKG eq "graphics" && $file =~ /[^m]text\.R$/) { next; }
    open(FILE, "< $file") or die "file $file cannot be opened";
    while (<FILE>) {
	$have_examples = 1
	    if ((/_ Examples _/o) || (/### \*+ Examples/));
	next if /^#/; # need to skip comment lines
	$have_par = 1 if (/[^a-zA-Z0-9.]par\(/o || /^par\(/o);
	$have_contrasts = 1 if /options\(contrasts/o;
    }
    close(FILE);
    if ($have_examples) {
	print "cleanEx(); nameEx(\"$nm\")\n";
    }

    print "### * $nm\n\n";
    print "flush(stderr()); flush(stdout())\n\n";
    open(FILE, "< $file") or die "file $file cannot be opened";
    my $dont_test = 0;	
    while (<FILE>) {
	## process \donttest
	$dont_test = 1 if /^## No test:/;
	print $_ unless $dont_test;
	$dont_test = 0 if /^## End\(No test\)/;
    }
    close(FILE);

    if($have_par) {
	## if there were 'par()' calls, now reset them:
	print "graphics::par(get(\"par.postscript\", pos = 'CheckExEnv'))\n";
    }
    if($have_contrasts) {
	## if contrasts were set, now reset them:
	print "options(contrasts = c(unordered = \"contr.treatment\"," .
	    "ordered = \"contr.poly\"))\n";
    }

}

### * Footer
print <<_EOF_;
### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\\\(> \\\\)?### [*]+" ***
### End: ***
quit('no')
_EOF_
