#/usr/bin/bash
mkdir temp
cd temp
rm -rf *.pdb
cp ../L/L_template.txt ref
id=$(<ref)
pdbid=${id:0:4}
echo $pdbid
#wget http://www.rcsb.org/pdb/files/$pdbid.pdb.gz
#gzip -d $pdbid.pdb.gz
wget http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/pdb/$pdbid/
mv index.html ref.pdb
cp ../L/*1.pdb L.pdb
cp ../H/*1.pdb H.pdb
cp ../TMalign .
./TMalign L.pdb ref.pdb -o LL
./TMalign H.pdb ref.pdb -o HH
cat HH_all_atm.pdb |awk '{print $0}'|awk '$5=="A"'   >TH
awk '{if ($5 == "A") gsub(/A/,"H",$5); printf "%4s%7.0f%5s%4s%2s%4.0f%12.3f%8.3f%8.3f%6.2f%6.2f\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}'<TH >THH
cat LL_all_atm.pdb |awk '{print $0}'|awk '$5=="A"'   >TL
awk '{if ($5 == "A") gsub(/A/,"L",$5); printf "%4s%7.0f%5s%4s%2s%4.0f%12.3f%8.3f%8.3f%6.2f%6.2f\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}'<TL >TLL
rm -rf candidate.pdb
cat THH >>candidate.pdb
cat ../TER >>candidate.pdb
cat TLL >>candidate.pdb
cat ../TER >>candidate.pdb
mkdir -p ../result
cp candidate.pdb ../result/
cd ..
rm -rf L
rm -rf H
rm -rf i*
rm -rf v
echo "job is down by FABbuilder01"
