#!perl
# Materials Studio's own atom types (and force-field charges) for every .mol2 in the current directory, with each force
# field named on the command line (Forcite's names: pcff, compass, cvff), the molecules read from the subfolder in/ (a file
# inside the script's own folder already belongs to the scripting project, and Import refuses it). Writes
# out/NAME.ms_FF.types: one line per atom in file order, "type charge". Run with RunMatScript from the folder above in/.
use strict;
use MaterialsScript qw(:all);
my $in = "D:/caps-build/msin";             # outside the scripting project folder (Import refuses files inside it)
if (@ARGV && $ARGV[0] =~ /^in=(.*)$/) { $in = $1; shift @ARGV; }
my @ffs = @ARGV ? @ARGV : ("pcff", "compass", "cvff");
use Cwd;
my $here = getcwd();
mkdir("$here/out");
opendir(my $dh, $in) or die "no $in";
my @files = sort grep { /\.mol2$/ && !/\.caps\./ } readdir($dh);
closedir($dh);
foreach my $f (@files) {
  (my $name = $f) =~ s/\.mol2$//;
  foreach my $ff (@ffs) {
    my $doc;
    eval {
      my $imp = Documents->Import("$in/$f");
      $doc = Documents->New("$name\_$ff.xsd");     # Forcite runs on an .xsd document only
      $doc->CopyFrom($imp);
      $imp->Discard;
      Modules->Forcite->Energy->Run($doc, Settings(CurrentForcefield => $ff, ChargeAssignment => "Forcefield assigned",
                                                    "3DPeriodicElectrostaticSummationMethod" => "Atom based"));
      open(my $o, ">", "$here/out/$name.ms_$ff.types") or die;
      foreach my $a (@{$doc->Atoms}) { printf $o "%s %.4f\n", $a->ForcefieldType, $a->Charge; }
      close($o);
    };
    print "$name $ff: ", ($@ ? "FAILED $@" : "ok"), "\n";
    eval { $doc->Discard if $doc; };
  }
}
print "done\n";
