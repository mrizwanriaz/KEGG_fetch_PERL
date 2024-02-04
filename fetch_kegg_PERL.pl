use strict;
use warnings;
use Carp;
use Cwd;
use File::Path qw(make_path);
use LWP::Simple;
use POSIX qw(strftime);

use Bio::KEGG::API;

# Create Bio::KEGG::API object
my $api = Bio::KEGG::API->new();

# Set up directory for saving output
my $dir = strftime("%Y%m%d", localtime(time)) . "//KEGG";
$dir = getcwd . '//' . $dir;
make_path($dir) unless -d $dir;

# List of KEGG databases
my @kegg_DBs = qw(pathway compound reaction rclass enzyme disease organism ko);

# Print start time
print "Start Time: " . strftime("%Y-%m-%d %H:%M:%S\n", localtime(time)) . "\n";

# Function calls
fetchKEGG_DBs(@kegg_DBs);
fetchPathwayLinkedEntries();
fetchReacionLinkedEntries();
fetchOrganismLinkedEntries();
fetchRPAIRs();

# Print end time and completion message
print "End Time: " . strftime("%Y-%m-%d %H:%M:%S\n", localtime(time)) . "\n";
print "All file(s) saved to $dir \n done... ";

# Sub-routines

# Fetch entries for specified KEGG databases
sub fetchKEGG_DBs {
    foreach (@_) {
        my @entries = $api->entry_list(database => $_);
        print scalar @entries . " entries fetched for $_ \n";
        writeFile("$_" . ".txt", @entries);
    }
}

# Fetch linked entries for pathways, including compounds, reactions, enzymes, and KO
sub fetchPathwayLinkedEntries {
    my $listFile = $dir . '//pathway.txt';
    open(my $refFile, '<:encoding(UTF-8)', $listFile) or die "Could not open file '$listFile' $!";

    my $cpdFile;
    open $cpdFile, '>>', $dir . '//pathwayCompounds.txt';
    my $rxnFile;
    open $rxnFile, '>>', $dir . '//pathwayReactions.txt';
    my $enzFile;
    open $enzFile, '>>', $dir . '//pathwayEnzymes.txt';
    my $koFile;
    open $koFile, '>>', $dir . '//pathwayKO.txt';

    my $kgmlDir = $dir . "/KGML_files";
    make_path($kgmlDir) unless -d $kgmlDir;

    while (my $row = <$refFile>) {
        chomp $row;
        my @columns = split /\t/, $row;
        my $id = substr($columns[0], 5, 12);
        my @array;

        # Fetch reactions
        @array = $api->linked_entries(target => 'rn', source => $id);
        print $rxnFile join("\r\n", @array);

        # Fetch compounds
        @array = $api->linked_entries(target => 'cpd', source => $id);
        print $cpdFile join("\r\n", @array);

        # Fetch enzymes
        @array = $api->linked_entries(target => 'enzyme', source => $id);
        print $enzFile join("\r\n", @array);

        # Fetch KO entries
        @array = $api->linked_entries(target => 'ko', source => $id);
        print $koFile join("\r\n", @array);

        # Fetch reference KGML
        open my $kgmlFile, '>>', "$kgmlDir/$id.kgml";
        my $kgml = fetchKGML(substr($id, 3, 8));
        print $kgmlFile $kgml if defined $kgml;
        close $kgmlFile;
        last;
    }

    close $cpdFile;
    close $rxnFile;
    close $enzFile;
    close $refFile;
}

# Fetch KGML file from URL because BIO::KEGG does not have any function for this
sub fetchKGML {
    my $url = 'http://rest.kegg.jp/get/rn' . $_[0] . '/kgml';
    return get($url);
}

# Fetch linked entries for organisms, including genes, enzymes, and KO
sub fetchOrganismLinkedEntries {
    my $listFile = $dir . '//organism.txt';
    open(my $refFile, '<:encoding(UTF-8)', $listFile) or die "Could not open file '$listFile' $!";

    my $geneFile;
    open $geneFile, '>>', $dir . '//organismGenes.txt';
    my $enzFile;
    open $enzFile, '>>', $dir . '//organismEnzymes.txt';
    my $koFile;
    open $koFile, '>>', $dir . '//organismKO.txt';

    while (my $row = <$refFile>) {
        chomp $row;
        my @columns = split /\t/, $row;
        my $id = $columns[1];
        my @array;

        # Fetch genes
        @array = $api->entry_list(database => $id);
        print $geneFile join("\r\n", @array);

        # Fetch enzymes
        @array = $api->linked_entries(target => 'enzyme', source => $id);
        print $enzFile join("\r\n", @array);

        # Fetch KO entries
        @array = $api->linked_entries(target => 'ko', source => $id);
        print $koFile join("\r\n", @array);
        last;
    }

    close $koFile;
    close $geneFile;
    close $enzFile;
    close $refFile;
}

# Fetch linked entries for reactions, including enzymes, and KO
sub fetchReacionLinkedEntries {
    my $enzFile;
    open $enzFile, '>>', $dir . '//reactionEnzymes.txt';
    my $koFile;
    open $koFile, '>>', $dir . '//reactionKO.txt';
    my $cpdFile;
    open $cpdFile, '>>', $dir . '//reactionCompound.txt';

    my @array;

    # Fetch compounds
    @array = $api->linked_entries(target => 'cpd', source => 'rn');
    print $cpdFile join("\r\n", @array);

    # Fetch enzymes
    @array = $api->linked_entries(target => 'enzyme', source => 'rn');
    print $enzFile join("\r\n", @array);

    # Fetch KO entries
    @array = $api->linked_entries(target => 'ko', source => 'rn');
    print $koFile join("\r\n", @array);

    close $koFile;
    close $cpdFile;
    close $enzFile;
}

# Fetch RPAIRs from RCLASS
sub fetchRPAIRs {
    my $listFile = $dir . '//rclass.txt';
    open(my $refFile, '<:encoding(UTF-8)', $listFile) or die "Could not open file '$listFile' $!";

    my $rpairFile;
    open $rpairFile, '>>', $dir . '//rpairs.txt';

    while (my $row = <$refFile>) {
        chomp $row;
        my @columns = split /\t/, $row;
        my $id = substr($columns[0], 3, 9);
        my $url = 'http://rest.kegg.jp/get/' . $id;

        my $res = get($url);
        if (defined $res) {
            print $rpairFile $res . "\n" . "###" . "\n";
        }
        last;
    }

    close $rpairFile;
    close $refFile;
}

# Utility function to write content to a file
sub writeFile {
    my $file = $dir . '//' . $_[0];

    if (-f $file) {
        unlink $file or croak "Cannot delete $file: $!";
    }

    open my $outFile, '>>', $file or croak "Cannot open $file: $OS_ERROR";

    foreach (@_[1 .. $#_]) {
        print $outFile "$_\n";
    }

    close $outFile or croak "Cannot close $file: $OS_ERROR";
}
