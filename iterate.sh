if [ $(ls -d [0-9][0-9][0-9][0-9] | tail -1) != "" ]
then
echo "directories exist"
lastdir=$(ls -d [0-9][0-9][0-9][0-9] | tail -1)
else
echo "no directories exist"
lastdir=-1
fi
echo $lastdir
lastdir=$(echo $lastdir | sed 's/^0*//')

newdir=$((++lastdir))
echo $newdir
mkdir $(printf "%04u" $newdir)
mkdir -p old
cp -r out/* old
mv data $(printf "%04u" $newdir)
mv in $(printf "%04u" $newdir)
mv out $(printf "%04u" $newdir)

