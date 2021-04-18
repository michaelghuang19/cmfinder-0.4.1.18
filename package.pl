use strict;
use Cwd;
# use File::Copy "cp";  # see note for 'cp' sub below
use File::Path qw(make_path);
use File::Basename;
use Getopt::Long;

my $initPwd=getcwd;
my $pdflatex="pdflatex";

my ($allowDirty,$overwriteDir);
my $packageRootDir="build-package";
my $infoPackageFileName="info.package";
my $packageVersion=undef;
my $verbose=0;
if (!GetOptions(
	 "allowDirty" => \$allowDirty, # DON'T USE THIS FOR RELEASES
	 "verbose" => \$verbose,
	 "overwriteDir" => \$overwriteDir
    )) {
    die "problem with command line params";
}

my $packageName=undef;
my $packageVersionVariable=undef;
my @autoreconfCommandList=();
my @noCopyFileList=();
my @symLinkFileList=(); # NOTE: in principle, we could scan for symbolic links and fix them, but I think this is a bit overly tricky, since (1) we need to calculate relative links, and (2) some links might be absolute (like 'depcomp').  So, I think it's easier to just do it more manually.
my @deleteAfterCopyList=();
my @pathChangeList=(); # [srcPath,destPath,configure.ac variable name (if any)]
my @editFileList=();
my @userManualFile=();
my @prefixLicenceFileList=();
my @additionalCopyIncludingSubdir=();
my %perlUnifyModuleSet=();
my %did_perlUnifyModuleSet=();
my $userManualUsesBibtex=1; # assume it does
my %perlRenameModuleSet=();

# read info.package, so we know what to do
if (!open(I,$infoPackageFileName)) {
    die "cannot open $infoPackageFileName";
}
while (<I>) {
    s/[\r\n]//g;
    s/\#.*$//; # deal with comments
    if (/^\s*$/) { # empty or whitespace-only lines are ignored
	next;
    }
    my $line=$_;
    my @f=split /\t/,$line;
    my $command=shift @f;
    my $okay=0;
    if ($command eq "packageName") {
	$okay=1;
	$packageName=shift @f;
    }
    if ($command eq "packageVersion") {
	$okay=1;
	$packageVersion=shift @f;
    }
    if ($command eq "packageVersionVariable") {
	$okay=1;
	$packageVersionVariable=shift @f;
    }
    if ($command eq "noCopyFileList") {
	$okay=1;
	push @noCopyFileList,@f;
    }
    if ($command eq "symLinkFileList") {
	$okay=1;
	my $linkDir=shift @f;
	my $targetDir=shift @f;
	my @l=@f;
	push @symLinkFileList,[$linkDir,$targetDir,\@l];
	#print "symLink: linkDir=$linkDir , targetDir=$targetDir , ".join(" , ",@l)."\n";
    }
    if ($command eq "pathChangeList") {
	$okay=1;
	push @pathChangeList,\@f;
    }
    if ($command eq "editFileList") {
	$okay=1;
	push @editFileList,@f;
    }
    if ($command eq "userManualFile") {
	$okay=1;
	push @userManualFile,\@f;
    }
    if ($command eq "autoreconfCommandList") {
	$okay=1;
	push @autoreconfCommandList,@f;
    }
    if ($command eq "userManualUsesBibtex") {
	$okay=1;
	$userManualUsesBibtex=shift @f;
	if ($userManualUsesBibtex ne "0" && $userManualUsesBibtex ne "1") {
	    die "in info.package line for userManualUsesBibtex: please use '0' or '1'";
	}
    }
    if ($command eq "prefixLicenceFile") {
	die "actually, I'm not implementing this.  I don't think it's really necessary, and it's just causing me to delay doing more important things";
	$okay=1;
	my $x={};
	$x->{file}=shift @f;
	@{$x->{patternList}}=@f;
	push @prefixLicenceFileList,$x;
    }
    if ($command eq "additionalCopyIncludingSubdir") { # for files that should be packaged, but that aren't stored in git
	$okay=1;
	my $subdir=shift @f;
	my $regex=shift @f;
	push @additionalCopyIncludingSubdir,{subdir=>$subdir,regex=>$regex};
    }
    if ($command eq "perlRenameModule") {
	$okay=1;
	my $perlFile=shift @f;
	my %renameMap=@f;
	%{$perlRenameModuleSet{$perlFile}}=%renameMap;
    }
    if ($command eq "perlUnifyModule") {
	$okay=1;
	my $perlFile=shift @f;
	my @moduleList=@f;
	$perlUnifyModuleSet{$perlFile}=\@moduleList;
	die "not implemented, sorry -- seemed like too much of a hassle -- I'll wait till people complain.  I went for the simpler perlRenameModule";
    }
    if (!$okay) {
	die "couldn't understand $command (line=$line)";
    }
}
close(I);

# check if required variables are defined
if (!defined($packageName)) {
    die "you must define the packageName in $infoPackageFileName";
}
if (!defined($packageVersion)) {
    die "you must define the packageVersion in $infoPackageFileName (set it to \"git\" to get the version from git)";
}
if (scalar @autoreconfCommandList==0) {
    die "you must define autoreconfCommandList in $infoPackageFileName (a likely default would be:\nautoreconfCommandList\tautoreconf -fvi";
}

# check if pathChangeList has a prefix of another, which screws up the logic
my $prevSrcPrefix=undef;
for my $spec (@pathChangeList) {
    my ($srcPrefix,$destPrefix,$varName)=@$spec;
    if (defined($prevSrcPrefix)) {
	if ($srcPrefix =~ /^$prevSrcPrefix/) {
	    die "in pathChangeLists, src dir \"$prevSrcPrefix\" is a prefix of \"$srcPrefix\".  Although it's perfectly reasonable to set things up like this, it's a hassle to implemented it, because the regular expressions in package.pl have difficulty with knowing what to match.  Please just add something to avoid the prefix problem.";
	}
    }
    $prevSrcPrefix=$srcPrefix;
}

for my $perlFile (keys %perlRenameModuleSet) {
    if (!EditFile($perlFile)) {
	die "file \"$perlFile\" was mentioned in perlRenameModule.  you must add this file to editFileList";
    }
}

# get the files from git
my @gitFileList=();
my $gitLsCmd="git ls-files|";
if (!open(G,$gitLsCmd)) {
    die "cannot open $gitLsCmd";
}
while (<G>) {
    s/[\r\n]//g;
    push @gitFileList,$_;
}
if (!close(G)) {
    die "problem reading from $gitLsCmd";
}

# even if the user hard-codes the version number, we still want to ask git if files are dirty, and we also want to record the git version in the destination
my $gitVersionCmd="git describe --dirty --always";
my $gitVersion=`$gitVersionCmd`;
if ($? != 0) {
    die "problem with $gitVersionCmd";
}
$gitVersion =~ s/[\r\n]//g;
if (($gitVersion =~ /dirty/) && !$allowDirty) {
    die "your current files are dirty (there are uncommitted changes).  Please commit changes or run this script with -allowDirty (for temporary testing only)";
}
my $gitHashCmd="git describe --abbrev=0 --always";
my $gitHash=`$gitHashCmd`;
if ($? != 0) {
    die "problem with $gitHashCmd";
}
$gitHash =~ s/[\r\n]//g;
my $gitDateCmd="git log -1 --format=\%cd --date=short";
my $gitDate=`$gitDateCmd`;
if ($? != 0) {
    die "problem with $gitDateCmd";
}
$gitDate =~ s/[\r\n]//g;

if ($packageVersion eq "git") {
    $packageVersion=$gitVersion;
}
if ($packageVersion eq "git-plusdate") {
    $packageVersion="$gitVersion-$gitDate";
}

my $baseName="$packageName-$packageVersion";
my $destDir="$packageRootDir/$baseName";
my $tarFileBase="$baseName.tgz";
my $tarFile="$packageRootDir/$tarFileBase";

if (-e $destDir) {
    if (!$overwriteDir) {
	die "destination dir \"$destDir\" exists, so I'm concerned that if files are supposed to be deleted, they would now be removed.  You can re-run with -overwriteDir to ignore this issue";
    }
}

my @srcFileList=@gitFileList;

for my $x (@additionalCopyIncludingSubdir) {
    my $cmd="find $x->{subdir} -regex \"$x->{regex}\" -print |";
    my @l=();
    if (!open(F,$cmd)) {
	die "cannot open $cmd";
    }
    while (<F>) {
	s/[\r\n]//g;
	push @l,$_;
    }
    if (!close(F)) {
	die "problem reading from $cmd";
    }
    if (scalar @l==0) {
	die "additionalCopyIncludingSubdir $x->{subdir} $x->{regex} : resulted in zero files.  That must be wrong.";
    }
    push @srcFileList,@l;
}

my %isPathMade=();
my %didEditFile=();
my %alreadyProcessedFile=(); # be robust to duplicate listings of files

for my $srcFile (@srcFileList) {
    
    if ($alreadyProcessedFile{$srcFile}) {
	next;
    }
    $alreadyProcessedFile{$srcFile}=1;
    
    if (SkipFile($srcFile)) {
	if ($verbose) {
	    print "SkipFile $srcFile\n";
	}
	next;
    }
    my $relDestFile=CalcDestFileName($srcFile);
    my $destFile="$destDir/$relDestFile";

    my $destSubDir=dirname($destFile);

    if ($verbose && 0) {
	print "make_path($destSubDir)\n";
    }

    if (!$isPathMade{$destSubDir}) {
	make_path($destSubDir);
	$isPathMade{$destSubDir}=1;
    }
    
    if (EditFile($srcFile)) {
	$didEditFile{$srcFile}=1;
	if ($verbose) {
	    print "EditFile $srcFile $destFile\n";
	}
	my %thisPerlRenameModuleSet=();
	if (defined($perlRenameModuleSet{$srcFile})) {
	    %thisPerlRenameModuleSet=%{$perlRenameModuleSet{$srcFile}};
	}
	my $isConfigureAc=0;
	if ($srcFile =~ /configure[.]ac$/) {
	    $isConfigureAc=1;
	}
	my $isC=0;
	if ($srcFile =~ /[.]c$/) {
	    $isC=1;
	}
	if ($srcFile =~ /[.]cpp$/) {
	    $isC=1;
	}
	my $isPerl=0;
	if ($srcFile =~ /[.]pl$/) {
	    $isPerl=1;
	}
	my $isAutomake=0;
	if ($srcFile =~ /[.]am$/) {
	    $isAutomake=1;
	}
	my $isLatex=0;
	if ($srcFile =~ /[.]tex$/) {
	    $isLatex=1;
	}
	my $madeSomeChange=0;
	if (!open(IN,$srcFile)) {
	    die "cannot open $srcFile";
	}
	if (!open(OUT,">$destFile")) {
	    die "cannot open $destFile";
	}
	while (<IN>) {
	    my $line=$_;
	    my $origLine=$line;
	    my $gotIt=0;
	    if ($isConfigureAc) {
		$gotIt=1;
		for my $spec (@pathChangeList) {
		    my ($srcPrefix,$destPrefix,$varName)=@$spec;
		    if (/^AC_CONFIG_FILES/) {
			$line =~ s/\[$srcPrefix/\[$destPrefix/g;
		    }
		    if (/^$varName/) {
			$line =~ s/\"$srcPrefix/\"$destPrefix/g;
		    }
		}
		$line =~ s/^$packageVersionVariable=.*$/$packageVersionVariable=$packageVersion/;
		$line =~ s/^AC_INIT.([^,]+),([^,]+),/AC_INIT(\1,$packageVersion,/;
	    }
	    if ($isC) {
		$gotIt=1;
		$line =~ s/^\#define +$packageVersionVariable +\".*\"/\#define $packageVersionVariable \"$packageVersion\"/;
	    }
	    if ($isAutomake) {
		$gotIt=1;
		$line =~ s/^$packageVersionVariable=.*$/$packageVersionVariable=$packageVersion/;
	    }
	    if ($isPerl) {
		$gotIt=1;
		$line =~ s/^my +\$$packageVersionVariable=\".*\";/my \$$packageVersionVariable=\"$packageVersion\";/;
		if ($line =~ /^use +(.*);$/) {
		    my $module=$1;
		    my $newModule=$thisPerlRenameModuleSet{$module};
		    if (defined($newModule)) {
			$line="use $newModule;\n";
		    }
		}
	    }
	    if ($isLatex) { # latex variable is hardcoded as myversion
		$gotIt=1;
		$line =~ s/\\newcommand[{]\\myversion[}].0.[{].*[}]$/\\newcommand{\\myversion}[0]{$packageVersion}/;
	    }
	    if (!$gotIt) {
		die "you forgot to implement a handler for files like $srcFile";
	    }
	    if ($line ne $origLine) {
		$madeSomeChange=1;
	    }
	    print OUT $line;
	}
	if (!close(IN)) {
	    die "problem closing $srcFile";
	}
	if (!close(OUT)) {
	    die "problem closing $destFile";
	}
	if (!$madeSomeChange) {
	    die "file $srcFile was not modified by Edits.  That's suspicious.  \$packageVersionVariable=\"$packageVersionVariable\"";
	}
	clone_permissions($srcFile,$destFile);
    }
    else {
	if ($verbose) {
	    print "cp $srcFile $destFile\n";
	}
	if (-l $srcFile) {
	    my $linkTo=readlink($srcFile);
	    if (!defined($linkTo)) {
		die "readlink($srcFile) failed $? $!";
	    }
	    if (substr($linkTo,0,1) eq "/") {
		die "the file $srcFile is not a relative link, so this is unlikely to work, and package.pl doesn't implement this case";
	    }
	    if (!symlink($linkTo,$destFile)) {
		die "symlink($linkTo,$destFile) failed: $? $!";
	    }
	}
	else {
	    #print "cp($srcFile,$destFile)\n";
	    if (!cp($srcFile,$destFile)) {
		die "problem copying $srcFile -> $destFile";
	    }
	}
    }
    if ($perlUnifyModuleSet{$srcFile}) {
	$did_perlUnifyModuleSet{$srcFile}=1;

	# edit the destFile, so that it's possible to compose editFileList and perlUnifyModule
	my @moduleList=@{$perlUnifyModuleSet{$srcFile}};
	
	die "not further implemented, sorry";
    }
}

my @missingEditFileList=();
for my $f (@editFileList) {
    if (!$didEditFile{$f}) {
	push @missingEditFileList,$f;
    }
}
if (scalar @missingEditFileList>0) {
    die "the following files were listed in editFileList, but we never encountered the file in copying.  Did you get the name right?  The list: ".join(" ",@missingEditFileList);
}

my $gitVersionFileName="$destDir/git-version";
if (!open(GIT,">$gitVersionFileName")) {
    die "cannot open $gitVersionFileName";
}
print GIT "$gitVersionCmd : $gitVersion\n";
print GIT "$gitHashCmd : $gitHash\n";
print GIT "$gitDateCmd : $gitDate\n";
close(GIT);

for my $umf (@userManualFile) {
    my ($userManualDir,$userManualFileBase)=@$umf;
    ChDir("$destDir/$userManualDir");
    RunCmd("$pdflatex $userManualFileBase");
    if ($userManualUsesBibtex) {
	RunCmd("bibtex $userManualFileBase");
    }
    RunCmd("$pdflatex $userManualFileBase");
    RunCmd("$pdflatex $userManualFileBase");
    RunCmd("rm -f *.aux *.bbl *.blg *.idx *.log *.out *.toc");
    ChDir("$initPwd");
    my $srcFile="$destDir/$userManualDir/$userManualFileBase.pdf";
    my $destFile="$destDir/$userManualFileBase.pdf";
    if (!cp($srcFile,$destFile)) {
	die "problem copying $srcFile -> $destFile";
    }
    unlink($srcFile);
}

ChDir($destDir);

for my $cmd (@autoreconfCommandList) {
    RunCmd($cmd);
}

# we don't need this cache dir
RunCmd("rm -rf autom4te.cache");

ChDir($initPwd);

# files (well just one file) needed for autoconf/automake, but that we don't need in the tar archive
for my $relDestBase (@deleteAfterCopyList) {
    my $destFile="$destDir/$relDestBase";
    unlink($destFile);
}

if (1) {
    for my $spec (@symLinkFileList) {
	my ($linkDir,$targetDir,$fileListRef)=@$spec;
	for my $file (@$fileListRef) {
	    my $linkFile="$destDir/$linkDir/$file";
	    my $numUpDirs=scalar split /\//,$linkDir;
	    my @upDirList=".." x $numUpDirs;
	    my $targetFile=join("/",@upDirList)."/"."$targetDir/$file";
	    unlink($linkFile); # remove any symbolic link that's now within the wrong link structure
	    if (!symlink($targetFile,$linkFile)) {
		die "problem with creating symbolic link $linkFile -> $targetFile (destDir=$destDir , linkDir=$linkDir , file=$file , numUpDirs=$numUpDirs)";
	    }
	}
    }
}

ChDir($packageRootDir);
RunCmd("tar czf $tarFileBase $baseName");

ChDir($initPwd);

print "\n\nPackaging compelete\n";
print " - your tar file is '$tarFile'\n";
print " - you can 'rm -r $packageRootDir' once you move the tar file\n";
print "\n";

sub EditFile {
    my ($srcFile)=@_;
    #print "EditFile(testing) , $srcFile , ".join(" , ",@editFileList)."\n";
    for my $ef (@editFileList) {
	if ($srcFile eq $ef) {
	    return 1;
	}
    }
    return 0;
}

sub SkipFile {
    my ($srcFile)=@_;
    for my $skipPrefix (@noCopyFileList) {
	if ($srcFile =~ /^$skipPrefix/) {
	    return 1;
	}
    }
    return 0;
}

sub CalcDestFileName {
    my ($srcFile)=@_;
    my $destFile=$srcFile;
    for my $spec (@pathChangeList) {
	my ($srcPrefix,$destPrefix,$varName)=@$spec;
	$destFile =~ s/^$srcPrefix/$destPrefix/g;
    }
    return $destFile;
}

sub RunCmd    {
    my ($cmd)=@_;
    print "Running $cmd\n";
    my $status=system "$cmd";
    if ($status!=0) {
	chdir($initPwd);
	die "failed with status $status (according to perl 'system' command)";
    }
}

sub cp {
    my ($src,$dest)=@_;
    RunCmd("cp \"$src\" \"$dest\""); # louise still has an older version of File::Copy that doesn't actually preserve permissions
    return 1;
}

sub clone_permissions {
    my ($src,$dest)=@_;
    my $perm = (stat($src))[2] & 07777;
    chmod($perm,$dest);
}

sub ChDir {
    my ($dir)=@_;
    if (!chdir($dir)) {
	die "chdir $dir failed";
    }
}
