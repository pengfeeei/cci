#include<fstream>
#include<iostream>
#include<time.h>
#include<sstream>
#include <string>
#include<cmath>
#include "memory.h"
#include <map>
#include <vector>
#include <algorithm>
using namespace std;


class Pengfei_DIY_Func
{
	public:
		double StringToDouble(string d_str);
		int StringToInt(string d_str);
		string DoubleToString(double d_num);
		string IntToString(int d_num);
		string RepStrInit(string d_str1,string d_str2,int d_num);
		string ReplaceAllTender(string d_str,string d_old_value,string d_new_value);
		string ReplaceAllRude(string d_str,string d_old_value,string d_new_value);
		string Trim(string d_str);
		string Included(string d_str,string d_data[],int d_num);
		string Split(string d_str,const char * d_split,string d_each[],int &d_num);
		double Var_Wei_jack(double d_fx[],double d_fx_all,int d_chunks,int d_each_num[],int d_marker_num);
		int Allele_adjust(string d_ori_snp,string d_ref_snp);
		string Allele_convert(string d_snp);

};

class Markers 
{ 	
	private: 
 		string index;
		double rna_exp; 
 		double cv; 
	public: 
		Markers(string index, double rna_exp,double cv)
		{
			this->index=index;
			this->rna_exp=rna_exp;
			this->cv=cv;
		}
		string get_index()
 		{ 
			return index;
		}
		double get_rna_exp()
		{ 
			return rna_exp;
		}
		double get_cv()
		{ 
			return cv;
		}

};


bool comparer(Markers m_a, Markers m_b)  
{ 
        return (m_a.get_rna_exp() < m_b.get_rna_exp());  
}
bool comparer2(double a, double b)  
{ 
        return (a < b);  
}


//*****************************************
//*****************************************
//*****************************************

int main(int argc,char *argv[])
{
	vector<Markers> vec;
	int i=0,beginmark=0,j=0,k=0,m=0,n,l,c,pop;
	ofstream outfile;
	ifstream infile;
	string str_line;
	string infilename,outfilename;
	string str_temp[50000];
	int locus_num=0;
	string str_swap;
	char * p;
	const char * split = "\t";
	map<string,string> map_rsid;
	int count=0;		
	Pengfei_DIY_Func myfunc;
	int gene_num;
	string index[50000];	
	string outstr[100];
	string cell_sample[70000];
        int cell_index[70000];
	double *gene_ex[70000];
        for(i=0;i<70000;i++)
        {
                gene_ex[i]=new double[30000];
        }

	double gene_mean[50000];
	double gene_var[50000];
	double cv[50000];
	double zscore[50000];
	double sum=0;
	string ident[30000];
	int ident_sub_num;
        int sample_num;
        int cell_num[10];
        int start[10];
        int end[10];
        string sample[1000];

	j=0;
        infilename="get_obj_data/list_obj.txt";
        infile.open(&infilename[0],ios::in);
        while(!infile.eof())
        {
                str_line="";
                getline(infile,str_line);
                if(str_line.length()<1)continue;
                j++;
                sample[j-1]=str_line;
        }
        infile.close();infile.clear();
        sample_num=j;
	cout<<"sample num "<<j<<endl;      

	outfilename="run_cci.sh";
	outfile.open(&outfilename[0],ios::out);
	for(i=0;i<sample_num;i++)
	{
		outfile<<"nohup bash run_cci.sh "<<sample[i]<<" &"<<endl;
	}
	outfile.close();outfile.clear();
}
//*****************************************
//*****************************************
//*****************************************

double Pengfei_DIY_Func::StringToDouble(string d_str)
{
	double d_num;
	stringstream d_sstr;
	d_sstr<<d_str;
	d_sstr>>d_num;
	return d_num;
}

int Pengfei_DIY_Func::StringToInt(string d_str)
{
	int d_num;
	stringstream d_sstr;
	d_sstr<<d_str;
	d_sstr>>d_num;
	return d_num;
}

string Pengfei_DIY_Func::DoubleToString(double d_num)
{
	string d_str;
	stringstream d_sstr;
	d_sstr<<d_num;
	d_sstr>>d_str;
	return d_str;	
}

string Pengfei_DIY_Func::IntToString(int d_num)
{
	string d_str;
	stringstream d_sstr;
	d_sstr<<d_num;
	d_sstr>>d_str;
	return d_str;	
}

string Pengfei_DIY_Func::RepStrInit(string d_str1,string d_str2,int d_num)
{
	for(int d_i=0;d_i<d_num;d_i++)
	{
		d_str1+=d_str2;
	}
	return d_str1;
}

string Pengfei_DIY_Func::ReplaceAllTender(string d_str,string d_old_value,string d_new_value)
{   
    for(string::size_type d_pos(0);d_pos!= string::npos;d_pos+= d_new_value.length())
	{
		if((d_pos=d_str.find(d_old_value,d_pos))!= string::npos)d_str.replace(d_pos,d_old_value.length(),d_new_value);
		else break;
	}
	return d_str;
}

string Pengfei_DIY_Func::ReplaceAllRude(string d_str,string d_old_value,string d_new_value)
{   
	while(true)
	{
		string::size_type d_pos(0);
		if((d_pos=d_str.find(d_old_value))!=string::npos)
			d_str.replace(d_pos,d_old_value.length(),d_new_value);
		else break;
	}
	return d_str;
}

string Pengfei_DIY_Func::Included(string d_str,string d_data[],int d_num)
{
	int d_i=0,d_j=0;
	string d_mark="0";
	for(d_i=0;d_i<d_num;d_i++)
	{
		if(d_str==d_data[d_i]){d_mark="1";break;}
	}
	return d_mark;

}

string Pengfei_DIY_Func::Trim(string d_str)
{
	int d_len=d_str.length();
	int d_i=0,d_j=0;
	string d_str2;
	for(d_i=0;d_i<d_len;d_i++)
	{
		if(d_str[d_i]!=' ')break;
	}
	for(d_j=d_len-1;d_j>=0;d_j--)
	{
		if(d_str[d_j]!=' ')break;
	}
	d_str2.assign(d_str,d_i,d_j-d_i+1);
	return d_str2;
}

string Pengfei_DIY_Func::Split(string d_str,const char * d_split,string d_each[],int &d_num)
{   
	char * d_p;
	int d_k=0;
    d_p = strtok(&d_str[0],d_split); 
	while(d_p!=NULL)
	{ 
		d_each[d_k++]=d_p;
		d_p = strtok(NULL,d_split); 
	}
	return d_str;
}


double Pengfei_DIY_Func::Var_Wei_jack(double d_fx[],double d_fx_all,int d_chunks,int d_each_num[],int d_marker_num)
{
	double d_sum1=0;
	double d_sum2=0;
	int d_i=0;
	double d_h=0;
	int d_chunks_valid=0;
	for(int d_i=0;d_i<d_chunks;d_i++)
	{		
		if(d_each_num[d_i]!=0)d_chunks_valid++;
	}
	for(int d_i=0;d_i<d_chunks;d_i++)
	{	
		if(d_each_num[d_i]==0)continue;
		d_sum2+=(1-(double)d_each_num[d_i]/d_marker_num)*d_fx[d_i];
	}
	for(int d_i=0;d_i<d_chunks;d_i++)
	{
		if(d_each_num[d_i]==0)continue;
		d_h=(double)d_marker_num/d_each_num[d_i];
		d_sum1+=pow(d_h*d_fx_all-(d_h-1)*d_fx[d_i]-d_chunks_valid*d_fx_all+d_sum2,2)/(d_h-1);
	}
	return pow(d_sum1/d_chunks_valid,0.5);
}

string Pengfei_DIY_Func::Allele_convert(string d_snp)
{
	if(d_snp=="A")return "T";
	else if(d_snp=="T")return "A";
	else if(d_snp=="G")return "C";
	else if(d_snp=="C")return "G";
	else return "N";
}
	
int Pengfei_DIY_Func::Allele_adjust(string d_ori_snp,string d_ref_snp)
{
	if(d_ori_snp[0]==d_ref_snp[0]&&d_ori_snp[1]==d_ref_snp[1])return 1;//ok
	else if(d_ori_snp[0]==d_ref_snp[1]&&d_ori_snp[1]==d_ref_snp[0])
	{
		if(d_ori_snp[0]=='A'&&d_ori_snp[1]=='T')return -9;//unsure, need to be deleted
		else if(d_ori_snp[0]=='T'&&d_ori_snp[1]=='A')return -9;
		else if(d_ori_snp[0]=='C'&&d_ori_snp[1]=='G')return -9;
		else if(d_ori_snp[0]=='G'&&d_ori_snp[1]=='C')return -9;
		else return 2;//be ref strand, but disordered 
	}	
	else if(d_ori_snp=="AG")
	{
		if(d_ref_snp=="TC")return -1;
		else if(d_ref_snp=="CT")return -2;
		else return -9;
	}
	else if(d_ori_snp=="AC")
	{
		if(d_ref_snp=="TG")return -1;
		else if(d_ref_snp=="GT")return -2;
		else return -9;
	}
	else if(d_ori_snp=="TC")
	{
		if(d_ref_snp=="AG")return -1;//not ref strand,change to opposite allele
		else if(d_ref_snp=="GA")return -2;//not ref strand, change to opposite allele, also disordered.
		else return -9;
	}	
	else if(d_ori_snp=="TG")
	{
		if(d_ref_snp=="AC")return -1;
		else if(d_ref_snp=="CA")return -2;
		else return -9;
	}
	else if(d_ori_snp=="CA")
	{
		if(d_ref_snp=="GT")return -1;
		else if(d_ref_snp=="TG")return -2;
		else return -9;
	}
	else if(d_ori_snp=="CT")
	{
		if(d_ref_snp=="GA")return -1;
		else if(d_ref_snp=="AG")return -2;
		else return -9;
	}
	else if(d_ori_snp=="GA")
	{
		if(d_ref_snp=="CT")return -1;
		else if(d_ref_snp=="TC")return -2;
		else return -9;
	}
	else if(d_ori_snp=="GT")
	{
		if(d_ref_snp=="CA")return -1;
		else if(d_ref_snp=="AC")return -2;
		else return -9;
	}
	else return -10;//absolutely wrong
}



