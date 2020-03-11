/*
 * Utility.cpp
 *
 *  Created on: Apr 15, 2016
 *      Author: jmauricio
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include "Utility.h"
#include "Structures.h"

using namespace std;

Utility::Utility() {
	// TODO Auto-generated constructor stub

}

Utility::~Utility() {
	// TODO Auto-generated destructor stub
}

double **Utility::createMatrixDouble(int _lines,int _rows){

	double **matrix;

    try
      {
    	matrix= new double *[_lines];
      }
    catch (bad_alloc& ba)
        {
         cerr << "Falta de memoria para matriz= new double *[linhas] bad_alloc : " << ba.what() << endl;
        }

    if(!matrix)
      {
        cout<<"Memória insuficiente para criar a matriz"<<endl;
        exit(1);
      }

    for(int i=0;i<_lines;i++)
     {
        try
          {
        	matrix[i]=new double [_rows];
          }
        catch (bad_alloc& ba)
            {
             cerr << "Falta de memoria para matrix[i]=new double [colunas] bad_alloc : " << ba.what() << endl;
            }

            if(!matrix)
             {
              cout<<"Memória insuficiente para criar a matriz"<<endl;
              exit(1);
             }
     }


  return matrix;
}
void Utility::printMatrix(double **matrix,int lines,int rows){
    cout<<"Aqui Matriz!!!"<<endl;

	    for(int i=0;i<lines;i++){
	       cout<<matrix[i][0];
	       //printf("%lg",matrix[i][0]);
	       //printf("%.21lf",matrix[i][0]);
	       for(int j=1;j<rows;j++){
	    	  cout<<" "<<matrix[i][j];
	    	  //printf(" %lg",matrix[i][j]);
	    	  //printf(" %.21lf",matrix[i][j]);

	       }
	       cout<<endl;
	       //printf("\n");
	     }

}

void Utility::printMatrixSINR(double **matrix,int lines,int rows){
    cout<<"Aqui Matriz!!!"<<endl;

	    for(int i=0;i<lines;i++){
	       cout<<std::setprecision(21)<<matrix[i][0]<<" ";
	       //printf("%.21lf  ",matrix[i][0]);
	       for(int j=1;j<rows;j++){
	    	  cout<<std::setprecision(21)<<" "<<matrix[i][j];
	          //printf("%.21lf  ", matrix[i][j]);
	       }
	       cout<<endl;
	     }

}
void Utility::convertTableToMw(double **SINR,double **SINR_Mw,int lines,int rows){

  double result, b;

  //***rcout<<"\n\nMatriz SINR_Mw com dados convertidos para Milliwatts"<<endl;
	for(int i=0;i<10;i++){
		for(int j=0;j<4;j++){
			if(SINR[i][j]!=0){

				b= SINR[i][j]/10.0;// dBm dividido por 10
				result= pow(10.0,b);//Converte de DBm para mW
				//***result= result*(pow(10,12));//Multiplica por 10^12 para se ter pelo menos três casas decimais quando for arredondar
				//***result= floor(result);//Arredonda o valor


				//result= result/(pow(10,12));//Divide por 10^12 para converter de volta
				//result= 10*log10(result);//Faz a conversão de mW para DBm


				SINR_Mw[i][j]=result;
				//***rprintf("%.21Lf ", SINR_Mw[i][j]);

			}
			else{
				 SINR_Mw[i][j]=0;
				 //***rcout<<"  -   ";
			    }
		}
		//***rcout<<endl;
	}

}
double Utility::convertDBMtoMW(double _value){
	 double result= 0.0, b;

	b= _value/10.0;// dBm dividido por 10
	result= pow(10,b);//Converte de DBm para mW

	//printf("%.21Lf ", _value);

 return result;
}
void Utility::printTableMwToDBm(double **SINR_Mw,int lines,int rows){

  double result;

	cout<<"\n\nMatriz com dados transformados de volta para DBm"<<endl;
	for(int i=0;i<10;i++){
		for(int j=0;j<4;j++){
			if(SINR_Mw[i][j]!=0){
				result= SINR_Mw[i][j];
				result= result/(pow(10,12));//Divide por 10^12 para converter de volta
				result= 10*log10(result);//Faz a conversão de mW para DBm

				printf("%.1f ", result);
			}
			else{
				 cout<<"  -   ";
			}
		}
		cout<<endl;
	}

}
void Utility::printGraph(vector< vector<Interference> > _graph){

	cout<<"\n\n       ====== Graph Interference ======\n\n"<<endl;

	for(int i=0;i<(int)_graph.size();i++){
		cout<<i<<": ";
		for(int j=0;j<(int)_graph[i].size();j++){
         cout<<"("<<_graph[i][j].distanceConnections<<","<<_graph[i][j].valueInterference<<") ";
         //printf("\n   -->(%d, %d, %.21Lf, %.21Lf)",_graph[i][j].idReceptor,_graph[i][j].idSender,_graph[i][j].distanceConnections, _graph[i][j].valueInterference);
		}
		cout<<endl;
	}

}

void Utility::printSetChannels(const vector<Channel> &setChannels){

	for(unsigned i=0;i<setChannels.size();i++){

		if(setChannels[i].bandwidth==0){
		 cout<<"B= 20";
		}
		else if(setChannels[i].bandwidth==1){
			  cout<<"B= 40";
		}
		else if(setChannels[i].bandwidth==2){
			cout<<"B= 80";
		}
		else{
			 cout<<"B= 160";
		}

		cout<<", F="<<setChannels[i].frequency<<" --->";
		for(unsigned j=0;j<setChannels[i].listConnections.size();j++){
			cout<<" "<<setChannels[i].listConnections[j];
		}
		cout<<endl;
	}

}

void Utility::printChannelsSlotTimeGroups(const vector<SetSlotTimes> &_timeSlotsGroups){

	cout<<"\n====== Lista de Conexoes ======"<<endl;

	for(unsigned i=0; i<_timeSlotsGroups.size(); i++){

		if(_timeSlotsGroups[i].listSpectrum[0].listChannels.size()>0 || _timeSlotsGroups[i].listSpectrum[1].listChannels.size()>0 || _timeSlotsGroups[i].listSpectrum[2].listChannels.size()>0){
			cout<<"\nTime[ "<<i<<" ]"<<endl;

			for(unsigned j=0; j<_timeSlotsGroups[i].listSpectrum.size(); j++){

				if(_timeSlotsGroups[i].listSpectrum[j].listChannels.size()>0){

					cout<<"Espectro= "<<j<<endl;

					for(unsigned s=0; s<_timeSlotsGroups[i].listSpectrum[j].listChannels.size(); s++){
						cout<<"---> Canal["<<s<<"](F="<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].frequency<<", B="<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth<<"): ";

						for(unsigned c=0; c<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections.size(); c++){

							cout<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[c]<<", ";

						}
						cout<<endl;
					}
					cout<<endl;

				}
			}
		}
	}


}

void Utility::printSlotTimeGroups(const vector<SetSlotTimes> &_timeSlotsGroups){

	cout<<"\n====== Lista de Canais ======"<<endl;

	for(unsigned i=0; i<_timeSlotsGroups.size(); i++){

		if(_timeSlotsGroups[i].listSpectrum[0].listChannels.size()>0 || _timeSlotsGroups[i].listSpectrum[1].listChannels.size()>0 || _timeSlotsGroups[i].listSpectrum[2].listChannels.size()>0){
			cout<<"\nTime[ "<<i<<" ]"<<endl;

			for(unsigned j=0; j<_timeSlotsGroups[i].listSpectrum.size(); j++){

				if(_timeSlotsGroups[i].listSpectrum[j].listChannels.size()>0){

					cout<<"Espectro= "<<j<<endl;
					cout<<"---> ";

					for(unsigned s=0; s<_timeSlotsGroups[i].listSpectrum[j].listChannels.size(); s++){
						cout<<"Canal["<<s<<"](F="<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].frequency<<", B="<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth<<"); ";
					}
					cout<<endl;

				}
			}
		}
	}


}

void Utility::printSlotTimeGroupsChannels(const vector<SetSlotTimes> &_timeSlotsGroups,int _slotTime,int _spectrum){

	cout<<"\n====== Canais por Tempo= "<<_slotTime<<" Espectro= "<<_spectrum<<" ======"<<endl;
	cout<<"--->";
	for(unsigned i=0; i<_timeSlotsGroups[_slotTime].listSpectrum[_spectrum].listChannels.size(); i++){
		cout<<"Canal["<<i<<"](F="<<_timeSlotsGroups[_slotTime].listSpectrum[_spectrum].listChannels[i].frequency
				<<", B="<<_timeSlotsGroups[_slotTime].listSpectrum[_spectrum].listChannels[i].bandwidth<<"); ";
	}
	cout<<endl;

}

void Utility::writeSlotTimeGroupsChannels(const vector<SetSlotTimes> &_timeSlotsGroups,string nameFile,double fo_Star){

 ofstream histogramFile;

 histogramFile.open(nameFile.c_str());

 if(histogramFile.is_open()){
	//cout<<"Arquivo histograma criado com sucesso "<<nameFile<<" !!!"<<endl;
 }
 else{
	  cout<<"Nao foi possivel abrir o arquivo do histograma: "<<nameFile<<endl;
	  histogramFile.clear();
	  histogramFile.close();
     }


 histogramFile<<"FO Star= "<<fo_Star<<"  \n";

	for(unsigned i=0; i<_timeSlotsGroups.size(); i++){

		if(_timeSlotsGroups[i].listSpectrum[0].listChannels.size()>0 || _timeSlotsGroups[i].listSpectrum[1].listChannels.size()>0 || _timeSlotsGroups[i].listSpectrum[2].listChannels.size()>0){

			histogramFile<<"Time[ "<<i<<" ]= "<<i<<"  \n";

			for(unsigned j=0; j<_timeSlotsGroups[i].listSpectrum.size(); j++){

				if(_timeSlotsGroups[i].listSpectrum[j].listChannels.size()>0){

					histogramFile<<"\nEspectro= "<<j<<"\n";

					for(unsigned s=0; s<_timeSlotsGroups[i].listSpectrum[j].listChannels.size(); s++){

						histogramFile<<"     Canal["<<s<<"](F="<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].frequency<<", B="<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].bandwidth<<"): ";

						for(unsigned c=0; c<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections.size(); c++){

							histogramFile<<_timeSlotsGroups[i].listSpectrum[j].listChannels[s].listConnections[c]<<"  ";

						}
						histogramFile<<"\n";
					}
					histogramFile<<"\n";

				}
			}
		}

		histogramFile<<"\n";
	}

	histogramFile.clear();
	histogramFile.close();
}

void Utility::writeHistogramTableBandwidthRate(const vector<SetSlotTimes> &_timeSlotsGroups,string nameFile,int **histogramTable,TimeTableBandwidthRate histogramTimeTable[],int _timeSlots){

 ofstream histogramFileTable;

 histogramFileTable.open(nameFile.c_str());

 if(histogramFileTable.is_open()){
	//cout<<"Arquivo histograma criado com sucesso "<<nameFile<<" !!!"<<endl;
 }
 else{
	  cout<<"Nao foi possivel abrir o arquivo do histograma: "<<nameFile<<endl;
	  histogramFileTable.clear();
	  histogramFileTable.close();
     }

  
 
  
  
  histogramFileTable<<"MCS      0 1 2 3 4 5 6 7 8 9 NT\n";
  for(int i=0; i<4; i++){

    if(i==0){
	 histogramFileTable<<"20MHz    ";
    }
	else if(i==1){
	        histogramFileTable<<"40MHz    ";
    }
	else if(i==2){
	       histogramFileTable<<"80MHz    ";
    }
    else{
    	  histogramFileTable<<"160MHz   ";
    }

    histogramFileTable<<histogramTable[i][0];

	//histogramFileTable<<"         "<<histogramTable[i][0];

	for(int j=1; j<11; j++){
	 histogramFileTable<<" "<<histogramTable[i][j];
	}

	histogramFileTable<<"\n";
  }

 histogramFileTable<<"\n\n\n";


  histogramFileTable<<"====== Histogramas por Tempo ======\n";

  for(int l=0; l<_timeSlots; l++){

	  histogramFileTable<<"Time_"<<l<<"\n";
	  histogramFileTable<<"MCS      0 1 2 3 4 5 6 7 8 9 NT\n";
	  for(int i=0; i<4; i++){

   		if(i==0){
	 	  histogramFileTable<<"20MHz    ";
    	}
		else if(i==1){
	          histogramFileTable<<"40MHz    ";
    	}
		else if(i==2){
	          histogramFileTable<<"80MHz    ";
    	}
    	else{
    	  	 histogramFileTable<<"160MHz   ";
   	 	}

    
	   histogramFileTable<<histogramTimeTable[l].timeTable[i][0];


		  //histogramFileTable<<"       "<<histogramTimeTable[l].timeTable[i][0];

		  for(int j=1; j<11; j++){
			  histogramFileTable<<" "<<histogramTimeTable[l].timeTable[i][j];
		  }

		  histogramFileTable<<"\n";
	  }

	  histogramFileTable<<"\n\n\n";

  }

  histogramFileTable.clear();
  histogramFileTable.close();
}

