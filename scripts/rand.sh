# Generate a random integer number between X (max) and Y (min)
X=120
Y=60
DIFF=$(($X-$Y+1))
N=1
for i in `seq $N`
do
    SEED=$RANDOM # A UNIX OS generated pseudorandom integer in the range 0 - 32767
    R=$(($(($SEED%$DIFF))+$Y))
    echo "Put system to sleep for $R seconds"
    #sleep $R
done