#!env bash
set -o pipefail

build_dir='build/bullseye'
rm -rf ${build_dir}
mkdir -p ${build_dir}
pushd ${build_dir}
wget http://www.cpan.org/authors/id/T/TO/TODDR/XML-Parser-2.46.tar.gz
tar -xf XML-Parser-2.46.tar.gz
pushd XML-Parser-2.46
perl Makefile.PL EXPATLIBPATH=$CONDA_PREFIX/lib EXPATINCPATH=$CONDA_PREFIX/include
make
make install
popd

cpanm Bio::DB::Fasta
cpanm Text::NSP
cpanm Array::IntSpan
cpanm MCE

echo $(pwd)
git clone https://github.com/mflamand/Bullseye bullseye
cp -r bullseye/Code/* $CONDA_PREFIX/bin/
chmod u+x $CONDA_PREFIX/bin/parseBAM.pl
chmod u+x $CONDA_PREFIX/bin/gtf2genepred.pl
chmod u+x $CONDA_PREFIX/bin/Find_edit_site.pl
chmod u+x $CONDA_PREFIX/bin/summarize_sites.pl

popd
