thread_num=4 # 定义最大线程数
ls -l |grep ^d | awk '{print $NF}' >oldTimes # 获取已经合并了的时刻
cd processor*
ls -l |grep ^d | awk '{print $NF}' >../newTimes # 获取分块的所有时刻
cd ..
reconList=`grep -wxvf oldTimes newTimes` # 获得尚未合并的时刻
i=0
for time in $reconList;do
    i=$(($i+1))
    reconstructPar -time $time &
    if [ "$i" -eq "$thread_num" ];then
        # 如果 i 达到线程上限了，则 wait，直到上面的命令都结束，然后重新对 i 计数
        wait
        i=0
    fi
done
