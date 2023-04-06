#!env bash
set -o pipefail

mkdir -p env_build
pushd env_build
wget http://www.cpan.org/authors/id/T/TO/TODDR/XML-Parser-2.46.tar.gz 
tar -xf XML-Parser-2.46.tar.gz
cd XML-Parser-2.46 
perl Makefile.PL EXPATLIBPATH=$CONDA_PREFIX/lib EXPATINCPATH=$CONDA_PREFIX/include
make
make install
popd

cpanm Bio::DB::Fasta
cpanm Text::NSP
cpanm Array::IntSpan
cpanm MCE
