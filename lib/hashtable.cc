#include "khmer.hh"
#include "hashtable.hh"
#include <iostream>

using namespace khmer;
using namespace std;

//
// filter_fasta_file: filter and trims a FASTA file into a new one
//

void Hashtable::filter_fasta_file(const std::string &inputfile,
                                  const std::string &outputfile, 
                                  int minLength, 
                                  int threshold)
{
   string line;
   ifstream infile(inputfile.c_str());
   ofstream outfile;
   outfile.open(outputfile.c_str());
   int isRead = 0;
   string name;
   string seq;

   if (infile.is_open())
   {
      while(!infile.eof())
      {
         getline(infile, line);
         if (line.length() == 0)
            break;

         if (isRead)
         {
            seq = line;

            if (get_min_count(seq) <= minLength) {
               outfile << ">" << name << endl;
               outfile << seq << endl;
            }

            /*
            int numPos = seq.length() - Hashtable::_ksize + 1;

            int readAbund[numPos];

            int start;
            int stop;
            
            for (int i = 0; i < numPos; i++)
            {
               string kmer= seq.substr(i, Hashtable::_ksize);
               readAbund[i] = Hashtable::get_min_count(kmer);            
            }

            start = 0;
            for (int i = 0; i < numPos; i++)
            {
               if (readAbund[i] >= threshold)
                  break;
               else
                  start++;
            }

            stop = numPos - 1;
            for (int i = (numPos-1); i >= 0; i--)
            {
               if (readAbund[i] >= threshold)
                  break;
               else
                  stop--;
            }

            if ((stop - start + Hashtable::_ksize) > minLength)
            {
               string mySeq = seq.substr(start,(stop-start)+Hashtable::_ksize);
               outfile << ">" << name << endl;
               outfile << mySeq << endl;
            }
            */

            name.clear();
            seq.clear();
         }
         else
         {
            name = line.substr(1, line.length()-1);
         }

         isRead = isRead? 0 : 1;
      }
   }
  
   infile.close();
   outfile.close();
}


//
// consume_fasta: consume a FASTA file of reads
//

void Hashtable::consume_fasta(const std::string &filename)
{
   string line;
   ifstream infile(filename.c_str());
   int isRead = 0;

   if (infile.is_open())
   {
     while (!infile.eof())
     {
       getline(infile, line);

       if (isRead)
         Hashtable::consume_string(line);
       
       isRead = isRead? 0 : 1;
     }
  }
}

//
// consume_string: run through every k-mer in the given string, & hash it.
//

void Hashtable::consume_string(const std::string &s)
{
  const char * sp = s.c_str();
  unsigned int length = s.length();

#if 0
  const unsigned int length = s.length() - _ksize + 1;
  for (unsigned int i = 0; i < length; i++) {
    count(&sp[i]);
  }
#else
  unsigned int mask = 0;
  for (unsigned int i = 0; i < _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

  unsigned long long int h = 0; 
  unsigned long long int r = 0;
  
  _hash(sp, _ksize, &h, &r);

  unsigned long long int bin = 0;

  if (h < r)
    bin = h % _tablesize;
  else
    bin = r % _tablesize;

  if (_counts[bin] != MAX_COUNT)
    _counts[bin]++;

  for (unsigned int i = _ksize; i < length; i++) {
    short int repr = twobit_repr(sp[i]);

    // left-shift the previous hash over
    h = h << 2;

    // 'or' in the current nt
    h |= twobit_repr(sp[i]);

    // mask off the 2 bits we shifted over.
    h &= mask;

    // now handle reverse complement
    r = r >> 2;
    r |= (twobit_comp(sp[i]) << (_ksize*2 - 2));

    if (h < r)
      bin = h % _tablesize;
    else
      bin = r % _tablesize;

    if (_counts[bin] != MAX_COUNT)
      _counts[bin]++;
  }

#endif // 0
}


BoundedCounterType Hashtable::get_min_count(const std::string &s)
{
  const unsigned int length = s.length();
  const char * sp = s.c_str();
  BoundedCounterType min_count, count;

  unsigned int mask = 0;
  for (unsigned int i = 0; i < (unsigned int) _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

  unsigned long long int h = 0;
  unsigned long long int r = 0;
  
  _hash(sp, _ksize, &h, &r);

  if (h < r)
    min_count = this->get_count(h);
  else
    min_count = this->get_count(r);  

  for (unsigned int i = _ksize; i < length; i++) {
    // left-shift the previous hash over
    h = h << 2;

    // 'or' in the current nt
    h |= twobit_repr(sp[i]);

    // mask off the 2 bits we shifted over.
    h &= mask;

    // now handle reverse complement
    r = r >> 2;
    r |= (twobit_comp(sp[i]) << (_ksize*2 - 2));

    if (h < r)
      count = this->get_count(h);
    else
      count = this->get_count(r);
    
    if (count < min_count) {
      min_count = count;
    }
  }
  return min_count;
}

BoundedCounterType Hashtable::get_max_count(const std::string &s)
{
  const unsigned int length = s.length();
  const char * sp = s.c_str();
  BoundedCounterType max_count, count;

  unsigned int mask = 0;
  for (unsigned int i = 0; i < (unsigned int) _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

  unsigned long long int h = 0;
  unsigned long long int r = 0;

  _hash(sp, _ksize, &h, &r);

  if (h < r)
    max_count = this->get_count(h);
  else
    max_count = this->get_count(r);

  for (unsigned int i = _ksize; i < length; i++) {
    // left-shift the previous hash over
    h = h << 2;

    // 'or' in the current nt
    h |= twobit_repr(sp[i]);

    // mask off the 2 bits we shifted over.
    h &= mask;

    // now handle reverse complement
    r = r >> 2;
    r |= (twobit_comp(sp[i]) << (_ksize*2-2));

    if (h < r)
      count = this->get_count(h);
    else
      count = this->get_count(r);    

    if (count > max_count) {
      max_count = count;
    }
  }
  return max_count;
}

