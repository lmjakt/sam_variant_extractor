#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h> // for atoi
#include <sys/mman.h>  // for mmap
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>  // these three to give me stat...
#include <stdio.h>  // this may not be needed for open
#include <stdlib.h>
#include <fcntl.h>
#include <pthread.h>
#include <omp.h>

// Compile with :
// g++ -O2 -fopenmp -lpthread -o readSam readSam.cpp

/////////////////////////////////////////////////////////
// THIS VERSION ASSUMES PAIRED READS FOR EVERYTHING    //
//                                                     //
// AND THAT ALL MAPPED POSITIONS FOR A GIVEN PAIR OF   //
// SEQUENCES                                           //
// WILL BE PRESENT NEXT TO EACH OTHER IN THE FILE      //
/////////////////////////////////////////////////////////

// A small simple program that takes a SAM file and some set of sequence files
// reads in the sam file, and then makes some structures summarising the coverage
// and the positions of mismatches.

// we have a fair amount of memory to play with, so will not bother to be careful
// with memory access

// We may want to consider trying openMP on this later on, but for now, let's
// just run as it is.

///////////// IF THIS WORKS AT ALL I SHOULD REBASE IT WITH 
///////////// THE FUNCTIONS DECLARED I A HEADER AND SO ON.

// a map for making reverse complementary sequences
// This ought to be defined somewhere else, possible as
// part of an object containing the map.
std::map<char, char> complement;
// but these have to be defined within a function

// a set of the allowable cigare codes.
std::map<char, unsigned int> cigarCodes;

//////// This is real UGLY, BAD, AND ERROR PRONE ////////
//////// Define a global variantsMutex, initialise in the main function
//////// and use when we need to insert a new map
pthread_mutex_t variants_mutex;

// A function for splitting a string on a special character
// copied from read_cuff.cpp (see the R_Stuff directory)
void split_string(std::vector<std::string>& words, const char *text, size_t length, char sep)
{
  unsigned int start = 0;
  unsigned int end = 0;
  while(start < length){
    while( end < (length-1) && text[end] != sep )
      ++end;
    unsigned int l = (text[end] == sep) ? end - start : 1 + end - start;
    words.push_back( std::string( text + start, l ) );
    ++end;
    start = end;
  }
}

bool find_fasta_id(std::string& line, std::string& id)
{
  if(line[0] != '>')
    return(false);
  size_t pos = line.find_first_of(" \t\n");
  if(pos == std::string::npos){
    id = line.substr(1);
    return(true);
  }
  id = line.substr(1, pos-1);
  return(true);
}

// something to keep the data of a sam alignment
// We will not try to store all the SAM alignments at the same time
// but will instead use this to pass the data between functions
// We could do something like store an index of the byte positions
// of the start of an alignment in the file. Then if we want to
// we could go back and look at interesting ones (say everything
// with variants. But for now, I think we'll just leave it like this).
struct sam_alignment {
  // a subset of the fields defined by the SAM standard
  std::string qName;
  unsigned int flag;
  std::string rName;
  unsigned int pos;
  unsigned char mapQuality;
  std::string cigar;
  std::string seq;
  std::vector<unsigned char> qual;
  
  // derived qualities
  bool forwardStrand; // set from the flag
  // the positions to which bases in seq match
  // obtained by parsing the cigar string
  std::vector<int> matchPos;
  
  // using query based positions we can directly
  // access the variants and the associated qualities 
  // from the sequence above.
  std::vector<unsigned int> variantPos; 
  
  // a constructor from a line of the file
  sam_alignment(const char* buffer, size_t& buf_offset, size_t buf_size){
    size_t buf_end = buf_offset;
    while(buf_end < buf_size - 1 && buffer[buf_end] != '\n')
      ++buf_end;
    size_t line_length = (buffer[buf_end] == '\n') ? buf_end - buf_offset : 1 + buf_end - buf_offset;
    std::vector<std::string> words;
    split_string(words, buffer + buf_offset, line_length, '\t');
    buf_offset = buf_end + 1;
    if(words.size() >= 11){
      qName = words[0];
      flag = (unsigned int)atoi(words[1].c_str());
      rName = words[2];
      pos = (unsigned int)atoi(words[3].c_str());
      mapQuality = (unsigned char)atoi(words[4].c_str());
      cigar = words[5];
      seq = words[9];
      qual.resize(seq.size());
      for(unsigned int i=0; i < seq.size() && i < words[10].size(); ++i)
	qual[i] = ((unsigned char)words[10][i]) - 33; // assume phred33
      forwardStrand = (bool)(flag & 16);
      // we could call functions to set the matchPos here, but
      // since those functions are happy if they take the struct as an argument
      // it's not really an idea to call them here. Making them member functions
      // would work, or rewriting them to take lower level arguments. Hmm.
      // in any case I don't feel like writing these as member functions
      // let's leave them elsewhere. 
    }
    // if not true, we will have a lot of empty strings. We can check for those,
    // but it's not very good way to do things.
  }
};
  
// al1 is always positioned to the left (smaller position)
struct sam_alignment_pair {
  sam_alignment al1;
  sam_alignment al2;
  bool good;
  
  int seperation;
  int coverage;

  sam_alignment_pair(sam_alignment& al1, sam_alignment& al2) :
    al1(al1), al2(al2) {
    good = false;
    // very cursory check.
    if(al1.qName.size() && al2.qName.size())
      good = true;
    
    ////// WARNING: these will be incorrect if we have insertions
    ////// or deletions in the code. To handle those correctly we need
    ////// to handle the cigar string. That is more conveniently handled
    ////// in the BAM format, but here we'll need to deal with the sam.
    seperation = al2.pos - (al1.pos + al1.seq.size());
    coverage = (al2.pos + al2.seq.size()) - al1.pos;
  }
};

// the key is the variant nucleotide (i.e. A,C,G or T)
// the vector of unsigned chars holds the quality calls of
// the individual variant calls.
typedef std::map<char, std::vector<unsigned char> > SNP;
struct vmap_mutex_pair {
  pthread_mutex_t m;
  std::map<unsigned int, SNP> v;
};

typedef std::map<std::string, vmap_mutex_pair > variantMap;


void init_complement_dictionary(){
  complement['A'] = 'T';
  complement['C'] = 'G';
  complement['G'] = 'C';
  complement['T'] = 'A';
  complement['N'] = 'N';
  
  complement['a'] = 'T';
  complement['c'] = 'G';
  complement['g'] = 'C';
  complement['t'] = 'A';
  complement['n'] = 'N';
  // and just in case there are gaps 
  // we do as below.
  complement['*'] = '*';
  complement['-'] = '-';
}

void init_cigar_dictionary(){
  // 1 => push_back pos, pos++
  cigarCodes['M'] = 1;
  cigarCodes['='] = 1;
  cigarCodes['X'] = 1;
  // 2 => pos++
  cigarCodes['D'] = 2;
  cigarCodes['N'] = 2;
  cigarCodes['P'] = 2;
  // 3 => push_back -1
  cigarCodes['I'] = 3;
  cigarCodes['S'] = 3;
  // 4 => push back -1, pos++
  cigarCodes['H'] = 5;
  // sum( M, I, S, =, X ) = length of seq
}

// returns a vector containing the positions (offset from 0) of
// the template sequence that aligns to the sequence.
// Insertions in the query sequence are given an offset of -1
// as they do not have a logical template alignment position.
void parseCigar(sam_alignment& al){
  // the cigar string is given as either \d+M or just M
  // indicating \d+ numbers or 1 M where M is the cigar code.
  // here we have used -1 to indicate a location where there is
  // no direct equivalent position in the genome. However, we
  // should consider adding a -2 or different code to indicate
  // an insertion in seq as opposed to sequence padding.
  std::string number = "";
  unsigned int n = 1;
  al.matchPos.reserve(al.seq.size() + 10);
  unsigned int pos=al.pos;
  for(unsigned int i = 0; i < al.cigar.size(); ++i){
    if(al.cigar[i] <= '9' && al.cigar[i] >= '0'){
      number.append(1, al.cigar[i]);
      continue;
    }
    if(!cigarCodes.count(al.cigar[i])){
      std::cerr << "Unvalid cigar code: " << al.cigar[i] << " in string: "
		<< al.cigar << "  from seq with id: " << al.qName 
		<< "and rName: " << al.rName << std::endl;
      exit(1);
    }
    char code = al.cigar[i];
    n = number.size() ? atoi(number.c_str()) : 1;
    number = "";
    // the following code is rather ugly, but we can modify it later
    unsigned int opCode = cigarCodes[ code ];
    //    std::cout << "parseCigar Cigar: " << al.cigar << "  number: " << number 
    //	      << "  n: " << n <<  " code: " <<  code 
    //      << "  opCode: " << opCode << std::endl;
    // we could use a switch statement here:
    if(opCode == 1){
      for(unsigned int j=0; j < n; ++j){
	al.matchPos.push_back(pos); ++pos;
      }
      continue;
    }
    if(opCode == 2){
      pos += n;
      continue;
    }
    if(opCode == 3){
      for(unsigned int j=0; j < n; ++j)
	al.matchPos.push_back(-1);
      continue;
    }
    if(opCode == 4){
      for(unsigned int j=0; j < n; ++j)
	al.matchPos.push_back(-1);
      pos += n;
      continue;
    }
  }
  //  for(unsigned int i=0; i < al.matchPos.size(); ++i){
  // std::cout << "  " << al.matchPos[i] ;
  //}
  //std::cout << std::endl;
  // the length of al.matchPos should be equal to the length of the sequence
  // if that is not the case, then we have made a mistake and should die here
  if(al.matchPos.size() != al.seq.size()){
    std::cerr << "al.matchPos.size != al.seq.size. Programmer error" << std::endl;
    exit(1);
  }
} 

// a function to read in the sequences in a fasta file
std::map<std::string, std::string> readFasta(std::string fname)
{
  std::map<std::string, std::string> sequences;
  std::ifstream in(fname.c_str());
  if(!in.is_open()){
    std::cerr << "Unable to open: " << fname << std::endl;
    exit(1);
  }
  std::string line;
  std::string id = "";
  bool id_found = false;
  while(getline(in, line)){
    // check for an a new id
    if(find_fasta_id(line, id)){
      sequences[id] = "";
      id_found = true;
      continue;
    }
    if(!id_found)
      continue;
    
    sequences[id] += line;
  }
  return(sequences);
}

struct d_mutex_pair {
  std::vector<int> d;
  pthread_mutex_t m;
};

/// We may need to move the call to pthread_mutex_init to somewhere else as
/// I believe we should not copy p_thread_mutexes. However, I'm not sure
/// if the return here copies the address. From test.cpp, it seems like
/// there is no copy when we assign the map from the return value, but not
/// sure if this is guaranteed.
std::map<std::string, d_mutex_pair > init_density_map(std::map<std::string, std::string>& seq)
{
  std::map<std::string, d_mutex_pair > counts;
  for(std::map<std::string, std::string>::iterator it=seq.begin();
      it != seq.end(); ++it){
    d_mutex_pair dmp;
    dmp.d = std::vector<int>(it->second.size(), 0);
    pthread_mutex_init( &dmp.m, NULL );
    counts.insert(make_pair((*it).first, dmp) );
  }
  // which is where we hope that this doesn't copy the result, but rather reassigns
  // a pointer (I believe this happens)
  return(counts);
}

/*
void increment_density_map(std::string seqId, std::map<std::string, d_map_pair >& dmap,
			   unsigned int beg, unsigned int range_length)
{
  std::map<std::string, std::vector<int> >::iterator it = dmap.find(seqId);
  if(it == dmap.end())
    return;
  if(beg + range_length > (*it).second.size()){
    std::cerr << "Range longer than target: " << beg
	      << " -> " << beg + range_length << "  > " 
	      << (*it).second.size() << "  " << (*it).first << std::endl;
    return;
  }
  unsigned int end = beg + range_length;
  while(beg < end){
    (*it).second[ beg ]++;
    ++beg;
  }
  return;
}
*/

// this function makes use of al.matchPos 
void increment_density_map(sam_alignment& al, std::map<std::string, d_mutex_pair >& dmap)
{
  std::map<std::string, d_mutex_pair >::iterator it = dmap.find(al.rName);
  if(it == dmap.end())
    return;
  //  std::cout << "1.1" << std::endl;

  if(!al.matchPos.size())
    parseCigar(al);
  //std::cout << "1.2" << std::endl;
  
  // here we should lock the mutex
  pthread_mutex_lock(&it->second.m);

  for(unsigned int i=0; i < al.matchPos.size(); ++i){
    //    std::cout << "1.3\t" << i << "\t" << al.matchPos[i] << std::endl;
    if(al.matchPos[i] > 0 && al.matchPos[i] <= it->second.d.size())
      it->second.d[ al.matchPos[i]-1 ]++;
  }
  pthread_mutex_unlock(&it->second.m);

// here we should unlock the mutex
}


// this would be a suitable place to use exceptions, but what the hell
// we need to work out what we want to return, or possibly what we want
// to write out here. 
sam_alignment_pair readPair(const char* sam, size_t& offset, size_t buf_size)
{

  sam_alignment al1(sam, offset, buf_size);
  sam_alignment al2(sam, offset, buf_size);
  
  if(al1.rName != "*")
    parseCigar(al1);
  if(al2.rName != "*")
    parseCigar(al2);

  // then let's get the relevant sequence from the genome
  // and go through and count mismatches. Hopefully we'll find some
  return(sam_alignment_pair(al1, al2));
}

std::string reverseComplement(std::string& s){
  std::string rc(s);
  //  rc.reserve(s.size());
  // using reserve means that we do not initialise the string
  std::map<char, char>::iterator it;
  for(unsigned int i=0; i < s.size(); ++i){
    it = complement.find( s[i] );
    if(it == complement.end()){
      std::cerr << "Unknown residue encountered: " << s[i] << std::endl
		<< "in seq: " << s << std::endl;
      exit(1);
    }
    rc[ rc.size() - (i+1) ] = it->second;
    //    rc.append(1, it->second );
  }
  //  std::cout << s << std::endl;
  // std::cout << rc << std::endl;
  return(rc);
}

void findVariants(std::map<std::string, std::string>& genome, sam_alignment& al)
{
  std::map<std::string, std::string>::iterator ref = genome.find(al.rName);
  if(ref == genome.end()){
    std::cerr << "No reference sequence with id: " << al.rName << std::endl;
    return;
  }
  // al.matchPos must be set for us to do this
  if(!al.matchPos.size()){
    // either set it or die? it should be set, but we should perhaps be
    // a little flexible. So lets complain and set it.
    std::cerr << "findVariants: al.matchPos not set. Setting it from here." << std::endl;
    parseCigar(al);
  }
  // We now need to find the beginning and end of the alignment positions
  // these we need to obtain from the matchPos vector. However, we cannot simply
  // take the front and back as these can contain -1 values (in the case
  // of soft clipped regions.
  unsigned int i=0;
  int beg=-1;
  int end=-1;
  while( (beg = al.matchPos[i]) == -1 && i < al.matchPos.size())
    ++i;
  // Not completely sure if the following works, but if not it should
  // crash pretty quickly.
  // It looks sort of elegeant, but I'm not a huge fan as it's not that
  // easy to read.
  i=al.matchPos.size();
  while( i > 0 && (end = al.matchPos[--i]) == -1 );

  if(beg == -1 || end == -1){
    std::cerr << "Unable to find matching positions on the sequences." << std::endl
	      << "This is a bug in the code which should be fixed. Aborting here." << std::endl;
    exit(1);
  }
  // Next we need to make sure that the positions are present in the ref sequence.
  // remember that in the SAM format the sequences are counted from 1 not 0
  if(beg <= 0 || end >= ref->second.size() || end <= 0){
    std::cerr << "countVariance: Match positions outside of reference sequence: "
	      << beg << " -> " << end << " outside of 1 -> " << ref->second.size() << std::endl;
    exit(1);
  }
  // Note that we have checked somewhere above that al.seq.size == al.matchPos.size
  // but since we are not using proper encapsulation this cannot be guarnteed here.
  //  std::string seq = al.forwardStrand ? al.seq : reverseComplement(al.seq);
  // The sequence in the sam file is already reverse complemented, so 
  // we should not do antything with it.
  for(unsigned int i=0; i < al.matchPos.size(); ++i){
    if(al.matchPos[i] <= 0)
      continue;

    if(al.seq[i] != ref->second[ al.matchPos[i]-1 ])
      al.variantPos.push_back(i);
  }
}

void countVariants(sam_alignment& al, variantMap& variants){
  if(!al.variantPos.size())
    return;
  std::vector<unsigned char> quals;
  std::vector<unsigned int> pos;
  std::vector<char> nucs;
  for(unsigned int i=0; i < al.variantPos.size(); ++i){
    unsigned int j = al.variantPos[i];
    quals.push_back( al.qual[j] );
    pos.push_back(al.matchPos[j] );
    nucs.push_back(al.seq[j]);
  }
  variantMap::iterator it = variants.find( al.rName );
  if(it == variants.end()){
    pthread_mutex_lock(&variants_mutex);
    vmap_mutex_pair vmap_pair;
    //    std::map<unsigned int, SNP> tmp;
    variants.insert(make_pair(al.rName, vmap_mutex_pair()));
    pthread_mutex_unlock(&variants_mutex);
    //    variants.insert(make_pair(al.rName, tmp));
    it = variants.find(al.rName);
    pthread_mutex_init(&it->second.m, NULL);
  }
  /// Whoah. This is real ugly. Really should make a more reasonable struct
  /// rather than trying to do this with typdeffing complex map structures.
  std::map<unsigned int, SNP>::iterator iit;
  SNP snp;
  // This locking thing is expensive... 
  pthread_mutex_lock(&it->second.m);
  for(unsigned int i=0; i  < quals.size(); ++i){
    iit = it->second.v.find( pos[i] );
    if(iit != it->second.v.end()){
      iit->second[ nucs[i] ].push_back( quals[i] );
    }else{
      it->second.v.insert(make_pair(pos[i], snp));
      it->second.v[ pos[i] ].insert( make_pair(nucs[i], std::vector<unsigned char>(1, quals[i]) ));
    }
  }
  //  std::cout << "variants.size: " << variants.size() << " " << it->first << " size : " << it->second.v.size() << std::endl;
  pthread_mutex_unlock(&it->second.m);
}

void writeVariant(std::ofstream& out, std::string chrome,
		  unsigned int pos, SNP& snp, std::map<std::string, d_mutex_pair >& densityMap,
		  std::map<std::string, std::string>& refSeq )
{
  // write
  // chrom position ref_nuc read_count V1;V2;V3 c1;c2;c3 qs1;qs2;qs3 q1,q1,q1;q2,;q3
  //
  // where: ref_nuc is the reference nucleotide
  //        read_count is the total number of reads (from densityMap)
  //        V1;V2, etc give the variant nucleotides
  //        c1;c2, etc give the variant counts
  //        qs1;qs2, etc give the quality sums
  //        q1,q2, give the qualities of the individual calls
  
  // we should also consider doing something like giving the individual read names,
  // read positions and byte offsets in the SAM file. But first check whether we are getting reasonable
  // numbers
  //  std::cout << "!!" << std::endl;
  std::vector<char> vars;
  std::vector<int> varCounts;
  std::vector<int> qualSums;
  std::vector<std::vector<unsigned char> > qualities;
  //std::cout << "??" << std::endl;
  out << chrome << "\t" << pos << "\t" << refSeq[chrome][pos-1] << "\t" << densityMap[chrome].d[pos-1];
  //std::cout << "*" << std::endl;
  for(SNP::iterator it=snp.begin(); it != snp.end(); ++it){
    vars.push_back(it->first);
    varCounts.push_back(it->second.size());
    qualSums.push_back(0);
    qualities.push_back(std::vector<unsigned char>());
    for(unsigned int i=0; i < it->second.size(); ++i){
      qualSums.back() += it->second[i];
      qualities.back().push_back( it->second[i] );
    }
  }
  //  std::cout << "@" << std::endl;
  out << "\t";
  for(unsigned int i=0; i < vars.size(); ++i)
    out << vars[i] << ";";
  out << "\t";
  for(unsigned int i=0; i < varCounts.size(); ++i)
    out << varCounts[i] << ";";
  out << "\t";
  for(unsigned int i=0; i < qualSums.size(); ++i)
    out << qualSums[i] << ";";
  out << "\t";
  for(unsigned int i=0; i < qualities.size(); ++i){
    for(unsigned int j=0; j < qualities[i].size(); ++j)
      out << (int)qualities[i][j] << ",";
    out << ";";
  }
  //  std::cout << "#" << std::endl;
  out << std::endl;
}

void writeDensityMap(std::string fname, std::map<std::string, d_mutex_pair >& dmap)
{
  std::ofstream out(fname.c_str());
  if(!out.is_open()){
    std::cerr << "Failed to open " << fname << " for writing: " << std::endl;
    return;
  }
  for(std::map<std::string, d_mutex_pair >::iterator it=dmap.begin();
      it != dmap.end(); ++it){
    unsigned int nameLength = it->first.size();
    unsigned int dataSize = it->second.d.size();
    out.write((const char*)&nameLength, sizeof(unsigned int));
    out.write((const char*)it->first.c_str(), it->first.size());
    out.write((const char*)&dataSize, sizeof(unsigned int));
    out.write((const char*)it->second.d.data(), sizeof(int) * it->second.d.size());
  }
  out.close();
}

//// That should be all the functions required.
//// In the main function do the following:

// 1. Read in the reference genome to a map<string, string>
// 2. Init the density map.
// 2. Open the SAM file, and read through until the last line of the
//    header. We are making an assumption here that we can identify
//    the last header line.
// 3. Read one pair at a time (readPair).
//    a. Increment the density map accordingly
//    b. detectVariants
//    
// Alignment pairs containing a variant are stored in an array of
// variants. Unfortunately we have no easy way of retrieving all the
// alignments for a given sequence. It might be possible to build
// some sort of index, (to byte positions in the file and then use
// fseek to the positions.
//
// given two indices, one for start positions and one for end positions
// might end up taking aroung 10GB of memory. So it's doable, but not
// very efficient. Well, as a multimap it might be a bit smaller, but,
// in any case lets simply count variant positions for the moment.
// make a map<chrom, map<pos, count> > to store the output positions.
// We can get the count of non-variant counts from the density map
// make sure to only count the actually observed bases for this!!

// WARNING: The following function does only very marginal error checking
//          It shouldn't fail, but if it does, I'm not really sure what the
//          program should do. And since it's unlikely we'll allow the program
//          to fail here.
std::string readNextId(const char* sam, size_t& offset, size_t max_offset)
{
  while(offset < max_offset && sam[offset] != '\n')
    ++offset;
  ++offset;  // to get to the beginning of the identifier.
  size_t id_end = offset;
  while(id_end < max_offset && sam[id_end] != '\t')
    ++id_end;
  while(id_end > offset && sam[id_end] != '_')
    --id_end;
  std::string id(sam + offset, id_end - offset);
  if(offset >= max_offset || id.size() < 2){
    std::cerr << "readNextId: fatal failure to read. Adjust proc_no and try again" << std::endl;
    exit(1);
  }
  return(id);
}

// inspect the logic of this carefully. It seems to be dangerous..
// Lets print out some stuff to see how it works..
void adjustOffset(const char* sam, size_t& offset, size_t max_offset){
  std::string id = readNextId(sam, offset, max_offset);
  // an empty while loop would seem to be sufficient to handle this
  while( readNextId(sam, offset, max_offset) == id );
}

// this function should be reentrant. Hence it needs to be called with
// an reference to its own variantMap. This is up to the caller to ensure.
// It is more difficult to ensure correct incrementation of the densityMap
// since to keep a single copy for each instance of each function is
// problematic. Similarly setting up mutexes to ensure stability seems
// excessively expensive. Is there a way around this?
// We can try to associate each sequence in the densityMap with a mutex
// and lock that when updating the map.

// maxMismatch is the maximum number of mismatches allowed for the sequence
// variants to be called.
void readSamSegment(const char* sam, size_t length, variantMap& variants, 
		    std::map<std::string, std::string>& refSeq,
		    std::map<std::string, d_mutex_pair>& densityMap,
		    int maxMismatch)
{
  size_t offset = 0;
  unsigned int alignCounter = 0;
  // only use information from alignments with the lowest number
  // of mismatches. This would be more elegant using a multimap
  // but that requires a value compare object, and I don't feel
  // like writing one at the moment
  std::map<unsigned int, std::vector<sam_alignment_pair> > alPairs;
  std::string lastId = "";
  while(offset < length){
    sam_alignment_pair al = readPair(sam, offset, length);
    if(!al.good){
      std::cerr << "Bad alignment " << std::endl;
      break;
    }
    if(al.al1.rName == "*" || al.al2.rName == "*")
      continue;
    // determine the seq id for both pairs
    size_t ss1 = al.al1.qName.find_last_of("_");
    size_t ss2 = al.al2.qName.find_last_of("_");
    if(ss1 == std::string::npos || ss2 == std::string::npos){
      std::cerr << "Unexpected id form: " << al.al1.qName << " or " << al.al2.qName << std::endl;
      exit(1);
    }
    std::string id1 = al.al1.qName.substr(0, ss1);
    std::string id2 = al.al2.qName.substr(0, ss2);
    if(id1 != id2){
      std::cerr << "Mismatch of pair ids: " << id1 << " and " << id2 << std::endl;
      exit(1);
    }
    alignCounter++;
    if(! (alignCounter % 100000) )
      std::cout << alignCounter << std::endl;
    // check for variants
    findVariants(refSeq, al.al1);
    findVariants(refSeq, al.al2);

    // if we haven't defined an id, or if the current id is still the same as the last
    // id then we go through here
    if(!lastId.size() || lastId == id1){
      alPairs[al.al1.variantPos.size() + al.al2.variantPos.size()].push_back(al);
      lastId = id1;
      continue;
    }
    
    // if we are here, lastId is defined, and is not equal to the current
    // id. Hence we need to first increment our counters.
    if(alPairs.size()){
      // increment the density map and countVariants for the
      // best scoring (smallest number of variants).
      std::map<unsigned int, std::vector<sam_alignment_pair> >::iterator it=alPairs.begin();
      // the following line works because C++ initialises the value to 0.
      // variantDistribution[it->first]++;
      if((int)it->first <= maxMismatch || maxMismatch < 0){
	for(std::vector<sam_alignment_pair>::iterator iit=it->second.begin();
	    iit != it->second.end(); ++iit){
	  increment_density_map(iit->al1, densityMap);
	  increment_density_map(iit->al2, densityMap);
	  
	  countVariants( iit->al1, variants);
	  countVariants( iit->al2, variants);
	}
      }
      alPairs.clear();
      alPairs.insert(make_pair( al.al1.variantPos.size() + al.al2.variantPos.size(),
				std::vector<sam_alignment_pair>(1, al) ));
      lastId = id1;
    }
  }
  
  // increment the last pair if it is defined.. 
  std::map<unsigned int, std::vector<sam_alignment_pair> >::iterator it=alPairs.begin();
  for(std::vector<sam_alignment_pair>::iterator iit=it->second.begin();
      iit != it->second.end(); ++iit){
    if((int)it->first <= maxMismatch || maxMismatch < 0){
      increment_density_map(iit->al1, densityMap);
      increment_density_map(iit->al2, densityMap);
      
      countVariants( iit->al1, variants);
      countVariants( iit->al2, variants);
    }
  }
 
  std::cout << "At the end of a readSamSegment" << std::endl;
  std::cout << "Variants size : " << variants.size() << std::endl;
}

// 
bool readSAM(const char *fname, std::map<std::string, std::string>& refSeq,
	     std::map<std::string, d_mutex_pair>& densityMap, variantMap& variants,
	     std::map<unsigned int, unsigned int>& variantDistribution, int maxMismatch)
{
  // Instead of opening and reading the file using file streams we're simply 
  // going to mmap. First make sure that it exists, and get the size of it
  // using stat.
  struct stat stat_buffer;
  int file_descriptor;
  size_t file_size;
  char* sam;
  size_t length;  // we get this using stat
  
  if( stat(fname, &stat_buffer) == -1 ){
    // no file ? 
    std::cerr << "Stat returned an error. Skipping " << fname << std::endl;
    return(false);
  }
  file_size = stat_buffer.st_size;

  if( (file_descriptor = open(fname, O_RDONLY)) == -1 ){
    std::cerr << "Failed to open " << fname << std::endl;
    return(false);
  }
  
  if(file_size == 0){
    std::cerr << "Obtained a file size of 0 for " << fname << std::endl;
    return(false);
  }
  
  // mmap the file
  std::cout << "Calling mmap" << std::endl;
  sam = (char*)mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, file_descriptor, 0);
  std::cout << "mmap returned" << std::endl;

  // we shouldn't need the file descriptor anymore so, lets' close it
  if(close(file_descriptor) == -1){
    std::cerr << "Unable to close file descriptor for " << fname << std::endl;
  }

  // first find the beginning of the reads. Find a line that does not begin
  // with @. Look for \n[^@]
  // we can be clever here to avoid going through every line, but it's probably not worth it..
  
  size_t align_start = 0;
  while(align_start < (file_size-1)){
    if(sam[align_start] == '\n' && sam[align_start + 1] != '@'){
      align_start++;
      break;
    }
    ++align_start;
  }
  /*
  std::cout << "The beginning of the alignments" << std::endl;
  fwrite(sam+align_start, 1, 300, stdout);
  std::cout << "\n" << std::endl; // forces a flush?
  */

  // That seems to work. We now wish to divide up the file into N chunks which will
  // be processed in parallel. Each chunk should have a beg and and end position.
  // we can store these in a vector. First make a simple estimte by dividing the 
  // area into equal chunks, then go to the begin and end points and make sure
  // to align them such that all entries for a given id are in the same location.

  unsigned int proc_no = 16;
  // we make a vector having proc_no + elements so that one vector
  // can contain both the begin and end positions.
  std::vector<size_t> begin_pos(proc_no + 1);  
  for(unsigned int i=0; i < proc_no; ++i)
    begin_pos[i] = align_start + (i * (file_size - align_start) / (proc_no + 1));

  begin_pos[proc_no] = file_size;

  // then refine the begin_pos
  
  for(unsigned int i=1; i < proc_no; ++i){
    adjustOffset(sam, begin_pos[i], begin_pos[i+1]);
  }

  /*
  std::cout << "fileSize:    " << file_size << std::endl;
  std::cout << "align_start: " << align_start << std::endl;
  for(unsigned int i=0; i < begin_pos.size(); ++i){
    std::cout << "\t" << i << "\t" << begin_pos[i] << std::endl;
    fwrite(sam + begin_pos[i]-300, 1, 300, stdout); std::cout << "\t|";
    fwrite(sam + begin_pos[i], 1, 35, stdout); std::cout << std::endl;
  }
  */
  //  std::vector<variantMap> vmaps(proc_no);
  // This function should be possible to parallelise,
  // but we have to make use of vmaps instead.. 
  //  omp_set_num_threads(6);
#pragma omp parallel for
  for(unsigned int i=0; i < proc_no; ++i){
    ///////    CHANGE variantMap to vmaps[i]  and do merging if we are going to
    ///////    parallelise this loop...
    std::cout << "calling readSamSegment" << std::endl;
    readSamSegment(sam + begin_pos[i], begin_pos[i+1] - begin_pos[i], variants,
		   refSeq, densityMap, maxMismatch);
  }

  // find some way to join the vmaps in parallel. This be too slow
  // and it might be better to have set of mutexes for the vmaps as well.
  // 


  // unmap the file
  if(munmap(sam, file_size) == -1){
    std::cerr << "Unable to unmap memory associated with " << fname << std::endl;
  }

  // this first for testing
  return(true);

}

int main(int argc, char **argv)
{
  init_complement_dictionary();
  init_cigar_dictionary();
  pthread_mutex_init( &variants_mutex, NULL );
  
  // we have some optional values
  int opt;
  int maxMismatch = -1;
  while( (opt = getopt(argc, argv, "m:")) != -1 ){
    switch(opt){
    case 'm':
      maxMismatch = atoi(optarg);
      break;
    default:
      std::cerr << "Unknown option : " << opt << std::endl;
      exit(1);
    }
  }
  std::cout << "optind is : " << optind << "\n"
	    << "argc      : " << argc << "\n"
	    << "maxMismatch is " << maxMismatch << "\n";
  if(optind > argc)
    std::cout << "argv[" << optind << "] = " << argv[optind] << std::endl;

  for(int i=optind; i < argc; ++i){
    std::cout << "\t" << i << " : " << argv[i] << std::endl;
  }
  // argument counting begins with the binary name
  if( (argc - optind) < 3){
    std::cerr << "usage readSam [-m maxMismatch] refFile.fa output_prefix alignments1.sam alignments2.sam ..." << std::endl;
    exit(1);
  }

  std::string filePrefix = argv[optind + 1];
  // the output_prefix shouldn't be file, if it is then probably the user did something wrong
  // we can use std::inline to a reply here..
  if(std::ifstream(filePrefix.c_str())){
    std::cerr << filePrefix << " specified as the output_prefix is a file that exists" << "\n"
	      << "This is probably not what you intended and to save you the\n"
	      << "trouble the program will abort here if you give the wrong answer\n"
	      << "type yes to continue: ";
    std::string reply;
    std::cin >> reply;
    if(reply != "yes")
      exit(1);
  }
    

  std::map<std::string, std::string> refSeq = readFasta(std::string(argv[optind]));

  // set up the density count:
  std::map<std::string, d_mutex_pair > densityMap = init_density_map(refSeq);

  // Then we have to decide how to store the variant information;
  // map< rName, map< pos, SNP> > 
  // SNP is a map< char, vector<qual> >
  // that would hold everything so that we can look at qualities of individual variations
  // the density.
  // we typedef this structure as a variantMap

  variantMap variants;
  std::map<unsigned int, unsigned int> variantDistribution;
  //  std::map<std::string, unsigned int> varCounts;

  // read all the sam files specified on the command line
  for(int i=(optind+2); i < argc; ++i){
    if(!readSAM(argv[i], refSeq, densityMap, variants, variantDistribution, maxMismatch))
      std::cerr << "Unable to read data from : " << argv[i] << std::endl;
  }

  // Then find some good way to write out the basic data.

  // this will be 4 times the size of the genome, plus a little bit
  std::cout << "Writing out density map" << std::endl;
  std::string densityMapFile = filePrefix + "_densityMap";
  writeDensityMap(densityMapFile, densityMap);
  
  std::cout << "Writing out variants map: " << variants.size() << std::endl;
  std::string variantsFile = filePrefix + "_variants";
  std::ofstream out(variantsFile.c_str());
  if(out.is_open()){
  // for each variant
    for(variantMap::iterator it = variants.begin(); it != variants.end(); ++it){
      //      std::cout << ".";
      for(std::map<unsigned int, SNP>::iterator iit=it->second.v.begin();
	  iit != it->second.v.end(); ++iit){
	//std::cout << "|" << std::endl;
	writeVariant(out, it->first, iit->first, iit->second, densityMap, refSeq);
	//std::cout << "|\n" << std::endl;
      }
    }
    //std::cout << std::endl;
  }else{
    std::cerr << "Unable to open file variants for writing" << std::endl;
  }
  out.close();
  // and let's write out the variant distribution;
  std::string vdFile = filePrefix + "_variantDistribution";
  std::ofstream vd(vdFile.c_str());
  if(vd.is_open()){
    for(std::map<unsigned int, unsigned int>::iterator it=variantDistribution.begin();
	it != variantDistribution.end(); ++it)
      vd << it->first << "\t" << it->second << "\n";
    vd.close();
  }
}

