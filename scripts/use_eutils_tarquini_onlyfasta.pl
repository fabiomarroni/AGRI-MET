#!/usr/local/bin/perl -w
use LWP::Protocol::https;
use LWP::Simple;

# Download protein records corresponding to a list of GI numbers.

$db = 'nucleotide';
$id_list = 'MH087439,MH087440,MH087441,MH087442,MH087443,MH087444,MH087445,MH087446,MH087447,KR528581.1,KM491305.1,KU194413.1,KX522755,KU312039.1,KF134124.1,KF134125.1,KF134123.1,KF686810.1,KT894101.1';

#assemble the epost URL
$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "epost.fcgi?db=$db&id=$id_list";

$output = get($url);

$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);

### include this code for ESearch-EFetch
#assemble the efetch URL
$url = $base . "efetch.fcgi?db=$db&query_key=$key&WebEnv=$web";
$url .= "&rettype=fasta&retmode=text";

#post the efetch URL
$data = get($url);
print "$data";