                            >  �     3  %   	  X    #   File:       tsLib.make#   Target:     tsLib#   Created:    Jago, Jul 2001, Stefano M. Iacus##   This script assumes that R binary already exists in #   folder 'src/macintosh/bin'##   Tested with:##   MPW Shell 3.6d7#   MrC C Compiler 5.0.0d3c1#   Universal Headers 3.4#   CarbonLib 1.4#DLib            = tsMAKEFILE        = {DLib}Lib.make�MondoBuild�    = {MAKEFILE}  # Make blank to avoid rebuilds when makefile is modifiedMacF2C     		=  "::::macintosh:Mac F2C"F2CDir          =  ::::macintosh:f2c:Includes        = -i ::::include: �				  -i ::::macintosh: �                  -i "{F2CDir}"				  Sym-PPC         = -sym offPPCCOptions     = {Includes} {Sym-PPC} -opt off -includes unix -w 35,2 -shared_lib_export on -d HAVE_CONFIG_H -d Macintosh �                  -d TARGET_API_MAC_CARBON=1 -prefix RHeaders.h -align power 				  ### Library directory ###LibDir     		= ::::macintosh:bin:library:{DLib}:libs### Source Files ### SrcFiles        =  	burg.c �					carray.c �					eureka.c �					filter.c �					init.c �					mburg.c �					myw.c �					pacf.c �					PPsum.c �					qr.c �					starma.c �					stl.cFortFiles		= 	eureka.f �					starma.f �					stl.f### Object Files ###ObjFiles-PPC    =   burg.o �					carray.o �					eureka.o �					filter.o �					init.o �					mburg.o �					myw.o �					pacf.o �					PPsum.o �					qr.o �					starma.o �					stl.o### Libraries ###LibFiles-PPC    =  �				  "{SharedLibraries}CarbonLib" �				  "{PPCLibraries}PPCCRuntime.o" �				  "{SharedLibraries}StdCLib" �				  "::::macintosh:bin:R"  ### Default Rules ### .o  �  .c  {�MondoBuild�}	{PPCC} {depDir}{default}.c -o {targDir}{default}.o {PPCCOptions}### Build Rules ###{DLib}Lib  ��  {ObjFiles-PPC} {LibFiles-PPC} {�MondoBuild�} #create export table    if `Exists :expvar`	 delete :expvar	end    catenate �.x > expvar# checks if modules directory exsists    if ! `Exists -d "{LibDir}"`	 echo "Creating libs directory" "{LibDir}"      NewFolder "{LibDir}"    end# Build the vfonts module		PPCLink �		-o {LibDir}:{DLib}Lib �		{ObjFiles-PPC} �		{LibFiles-PPC} �		{Sym-PPC} �		-mf -d �		-t 'shlb' �		-c '????' �		-xm s �		-@export expvar# Removing mass    Delete   :�.x	Delete   :�.o     Delete   :expvar	### Required Dependencies ### dblcen.c  � dblcen.odistance.c  � distance.ohclust-utils.c � hclust-utils.o   hclust.c  � hclust.oinit.c  � init.okmns.c  � kmns.o# obj dependenciesburg.o 		�	burg.ccarray.o	�	carray.ceureka.o	�	eureka.cfilter.o	�	filter.cinti.o		�	init.cmburg.o		�	mburg.cmyw.o		�	myw.cpacf.o		�	pacf.cPPsum.o		�	PPsum.cqr.o		�	qr.cstarma.o	�	starma.cstl.o		�	stl.c# f2c dependencieseureka.c � eureka.f	 {MacF2C} -wait {FortFiles} starma.c � starma.fstl.c � starma.f### Optional Dependencies ###### Build this target to generate "include file" dependencies. ###Dependencies  �  $OutOfDate	MakeDepend �		-append {MAKEFILE} �		-ignore "{CIncludes}" �		-objdir ":" �		-objext .o �		{Includes} �		{SrcFiles}     �   �   a                                                                                                                                                                                                                                                   J 	Monaco    5�@                     < $� < $�����  �  �  � ��    < $� < $�  �             P��)?�20         �������                iacus   
tsLib.make    mpw makefile      �   �   aSORT� �  R MPSR  ckid   *���        ���   N     �     p��Projector DataTEXTMPS  ����                  