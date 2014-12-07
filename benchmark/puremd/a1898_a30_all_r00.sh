#PBS -S /bin/bash
#PBS -l nodes=1,walltime=5000:00:00
 
cd /home/chengtao/bboost/1898/1898k/a30_all/r00
echo $HOSTNAME > host.log
mkdir /state/partition1/chengtao
mkdir /state/partition1/chengtao/a1898_a30_all_r00
cp geo /state/partition1/chengtao/a1898_a30_all_r00/
cp ffield /state/partition1/chengtao/a1898_a30_all_r00/
cp ffield.ext /state/partition1/chengtao/a1898_a30_all_r00/
cp control /state/partition1/chengtao/a1898_a30_all_r00/
cd /state/partition1/chengtao/a1898_a30_all_r00/
/home/chengtao/bboost/code/SerialReax geo ffield control 
gzip *
 
cp * /home/chengtao/bboost/1898/1898k/a30_all/r00
if [ $status != 0 ]
then
        echo "Copying results from `hostname`:/state/partition1/chengtao/a1898_a30_all_r00/ back to $PBS_O_WORKDIR failed." 
        echo "After fixing the problem be sure to remove the directory `hostname`:/state/partition1/chengtao/a1898_a30_all_r00/"
else
        cd  /home/chengtao/bboost/1898/1898k/a30_all/r00
        rm -rf /state/partition1/chengtao/a1898_a30_all_r00/
fi
 
