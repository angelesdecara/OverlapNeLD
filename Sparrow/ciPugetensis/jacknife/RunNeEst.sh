#!/bin/bash
#cp ~/Downloads/Zip\ Folder_64_bit_191125/Ne2-1L .
./Ne2-1L c:inputNeEstpop0.txt 
for i in {1..101};do
    j=$((i-1))
    sed -i "s/pop${j}/pop${i}/g" inputNeEstpop0.txt
    ./Ne2-1L c:inputNeEstpop0.txt
done
