echo "compiling ... "
g++ script/s1_dis.cpp -o bin/s1_dis
echo "g++ script/s1_dis.cpp -o bin/s1_dis"
g++ script/s2_dis_bf.cpp -o bin/s2_dis_bf
echo "g++ script/s2_dis_bf.cpp -o bin/s2_dis_bf"
g++ script/s3_cor.cpp -o bin/s3_cor
echo "g++ script/s3_cor.cpp -o bin/s3_cor"
g++ script/s4_cor_bf.cpp -o bin/s4_cor_bf
echo "g++ script/s4_cor_bf.cpp -o bin/s4_cor_bf"
g++ script/s6_embedding_dis.cpp -o bin/s6_embedding_dis
echo "g++ script/s6_embedding_dis.cpp -o bin/s6_embedding_dis"
g++ script/s7_embedding_dis_bf.cpp -o bin/s7_embedding_dis_bf
echo "g++ script/s7_embedding_dis_bf.cpp -o bin/s7_embedding_dis_bf"
g++ script/s8_combine_bf.cpp -o bin/s8_combine_bf
echo "g++ script/s8_combine_bf.cpp -o bin/s8_combine_bf"

# an example for visualization 
g++ script/s9_readscore_for_plot.cpp -o bin/s9_readscore_for_plot
echo "g++ script/s9_readscore_for_plot.cpp -o bin/s9_readscore_for_plot"

echo "done."



