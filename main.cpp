#define N 1024
#define NoError 0
#define OpenFileError 1
#define FileReadError 2
#define RandgeError 3

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <iostream>

#include "randomc.h"

int32_t error;


size_t L;//размер решётки
double T;//температура
double a;//параметр для задания начального состояния
size_t mcs_max;//максимальное время в шагах Монте-Карло
size_t conf_max;//максимальное количество конфигураций (экспериментов/опытов/прогонок)

int64_t seed;//зерно для генератора случайных чисел

int8_t sp[N][N];//двумерная решетка спинов

double *Magnetic;//намагниченность
double *Magnetic2;//квадрат намагниченности

double *_Magnetic_;//намагниченность
double *_Magnetic2_;//квадрат намагниченности

//int8_t размер 8 бит(1байт) аналогичен char, от 0 до 255 безнаковый, знаковый от -128 до 127. 0xff
//int16_t
//int32_t
//int64_t


int32_t Initial(void);

inline int8_t sosedu(int8_t inSp[N][N], size_t inI, size_t inJ);
inline int8_t soseduDiagonal(int8_t inSp[N][N], size_t inI, size_t inJ);

inline double w(const double inDE, const double inT);
inline double calculateMagnetic(int8_t inSp[N][N]);
inline double calculateEnergy(int8_t inSp[N][N]);
int main()
{
    char fname[256];
    FILE *outMagnetic, *outGeneral;
    double dE;//E1,E2,dE;
    size_t mcs,conf;//временные шаги и прогонки.
    double r;
    int rI,rJ;
    CRandomMersenne Mersenne(time(0));

    error=Initial();


    #if WIN64
        if(error){printf_s("\n Initialization error %d",error); return error;}
    #endif

    #if __linux
        if(error){printf("\n Initialization error %d",error); return error;}
    #endif

    seed=time(0)+(3*(mcs_max%Mersenne.IRandomX(1,11))*Mersenne.IRandomX(1,200))*Mersenne.IRandom(L,L*N);
    Mersenne.RandomInit(seed);

    Magnetic=(double *)malloc(mcs_max*sizeof(double));//выделяем памят под намагниченность
    Magnetic2=(double *)malloc(mcs_max*sizeof(double));

    _Magnetic_=(double *)malloc(mcs_max*sizeof(double));//выделяем памят под намагниченность
    _Magnetic2_=(double *)malloc(mcs_max*sizeof(double));


    for(mcs=0;mcs<mcs_max;++mcs)
    {
        _Magnetic_[mcs]=Magnetic[mcs]=0;
        _Magnetic2_[mcs]=Magnetic2[mcs]=0;
    }
    for(conf=0;conf<conf_max;++conf)
    {
        //задаем начальное состояние для каждой новой конфигурации
        for(size_t i=0;i<L;++i)
        {
            for(size_t j=0;j<L;++j)
            {
                r=Mersenne.Random();
                if(a>r)//r случайное число от 0 до 1. a -- параметр от 0 до 1. При этом а=1 все спины верх, а=0 все спины вниз. а=1/2 половина спинов верх половина спинов вниз
                {
                    sp[i][j]=1;
                }
                else
                {
                    sp[i][j]=-1;
                }
            }
        }

        //начинаем считать перевороты спинов от шагов монте-карло. один шаг МОНТЕ-КАРЛО РАВЕН ВЕРОЯТНОСТИ ПЕРЕВОРОТА ВСЕХ СПИНОВ!!!
        for(mcs=0;mcs<mcs_max;++mcs)
        {
            if(mcs%(mcs_max/10)==0)//выводим через каждые mcs_max/10 шагов время, чтобы убедится, что программа не зависла и работает.
            {
                #if WIN64
                    printf_s("\nconf=%lld\tmcs=%lld",conf,mcs);
                #endif

                #if __linux
                    printf("\nconf=%ld\tmcs=%ld",conf,mcs);
                #endif

//                if(mcs>0)
//                {
//                    double E1,E2,dE1,dE2;
//                    size_t ii=112;
//                    size_t jj=89;
//                    E1=calculateEnergy(sp);

//                    dE1=2*(sp[ii][jj]*sosedu(sp,ii,jj)+0.5*sp[ii][jj]*soseduDiagonal(sp,ii,jj));
//                    sp[ii][jj]=-sp[ii][jj];
//                    E2=calculateEnergy(sp);
//                    dE2=E2-E1;

//                    printf("\n dE1=%lf\t dE2=%lf",dE1,dE2);
//                    return 0;
//                }
            }
            //считываем намагниченность
            Magnetic[mcs]=fabs(calculateMagnetic(sp));
            Magnetic2[mcs]=Magnetic[mcs]*Magnetic[mcs];

            _Magnetic_[mcs]+=Magnetic[mcs];
            _Magnetic2_[mcs]+=Magnetic2[mcs];

            //пытаемся перевернуть случайный спин в решетке
            for(size_t l=0;l<L*L;++l)
            {
                rI=Mersenne.IRandomX(0,L-1);
                rJ=Mersenne.IRandomX(0,L-1);

                dE=2*sp[rI][rJ]*sosedu(sp,rI,rJ);
                //алгоритм тепловой бани
                r=Mersenne.Random();
                if(r<w(dE,T))
                {
                    sp[rI][rJ]=-sp[rI][rJ];
                }

            }
        }
        #if WIN64
            sprintf(fname,"outMagnetic_conf=%lld.dat",conf);
        #endif

        #if __linux
            sprintf(fname,"outMagnetic_conf//outMagnetic_conf=%ld.dat",conf);
        #endif

        outMagnetic=fopen(fname,"w+");
        if(!outMagnetic){printf("File \"outMagnetic\" no open");return OpenFileError;}
        for(mcs=0;mcs<mcs_max;++mcs)
        {
            #if WIN64
                fprintf_s(outMagnetic,"%lld\t%lf\t%lf\n",mcs,Magnetic[mcs]/conf_max,Magnetic2[mcs]/conf_max);
            #endif

            #if __linux
                fprintf(outMagnetic,"%ld\t%lf\t%lf\n",mcs,Magnetic[mcs],Magnetic2[mcs]);
            #endif

        }
        fclose(outMagnetic);
    }//закончили цикл по прогонкам

    #if WIN64
        sprintf(fname,"outMagnetic_T=%lf.dat",T);
    #endif

    #if __linux
        sprintf(fname,"outMagnetic_T=%lf.dat",T);
    #endif

    outMagnetic=fopen(fname,"w+");
    if(!outMagnetic){printf("File \"outMagnetic\" no open");return OpenFileError;}
    size_t count=0;
    double m,m2;
    m=m2=0;
    for(mcs=0;mcs<mcs_max;++mcs)
    {
        double wm,wm2;
        wm=_Magnetic_[mcs]/conf_max;
        wm2=_Magnetic2_[mcs]/conf_max;
        #if WIN64
            fprintf_s(outMagnetic,"%lld\t%lf\t%lf\n",mcs,wm,wm2);
        #endif

        #if __linux
            fprintf(outMagnetic,"%ld\t%lf\t%lf\n",mcs,wm,wm2);
        #endif

        if(mcs>mcs_max*0.8)
        {
            count++;
            m+=wm;
            m2+=wm2;
        }

    }
    fclose(outMagnetic);

    m=m/count;
    m2=m2/count;


    #if WIN64
        sprintf(fname,"outT_L=%lld.dat",L);
    #endif

    #if __linux
        sprintf(fname,"outT_L=%ld.dat",L);
    #endif

    outGeneral=fopen(fname,"a");
    if(!outGeneral){printf("File \"outGeneral\" no open");return OpenFileError;}

    fprintf(outGeneral,"%lf\t%lf\t%lf\t%lf\n",T,m,m2,L*L*(m2-m*m)/T);
    fclose(outGeneral);

    free(Magnetic);
    free(Magnetic2);
    free(_Magnetic_);
    free(_Magnetic2_);

    #if WIN64
        printf_s("\nThe End!\n");
    #endif

    #if __linux
        printf("\nThe End!\n");
    #endif


    return NoError;
}


#if WIN64
int32_t Initial(void)
{
    FILE *inFile;

    inFile=fopen("Config.dat","r");//открыли файл для считывания
    if(!inFile){printf_s("File \"Config.dat\" no open!"); return OpenFileError;}

    //считываем данные из файла и проверяем -- получилось ли считать данные.

    //считываем температуру
    if(fscanf_s(inFile,"T=%lf;\n",&T))
    {
        printf_s("T=%lf;\n",T);
    }
    else
    {
        printf_s("Data \"T\" from the file have not been read");
        return FileReadError;
    }

    //считываем параметр, который отвечает за начальное состояние
    if(fscanf_s(inFile,"a=%lf;\n",&a))
    {
        printf_s("a=%lf;\n",a);
    }
    else
    {
        printf_s("Data \"a\" from the file have not been read");
        return FileReadError;
    }

    //считываем размер системы
    if(fscanf_s(inFile,"L=%lld;\n",&L))//%lld для linux заменить на ld или d?
    {
        printf_s("L=%lld;\n",L);
        if(L>256)
        {
            printf_s("L should be no more than 256");
            return RandgeError;
        }
    }
    else
    {
        printf_s("Data \"L\" from the file have not been read");
        return FileReadError;
    }

    //считываем максимальное количество временных шагов
    if(fscanf_s(inFile,"mcs_max=%lld;\n",&mcs_max))
    {
        printf_s("mcs_max=%lld;\n",mcs_max);
    }
    else
    {
        printf_s("Data \"mcs_max\" from the file have not been read");
        return FileReadError;
    }

    //считываем максимальное количество конфигураций
    if(fscanf_s(inFile,"conf_max=%lld;\n",&conf_max))
    {
        printf_s("conf_max=%lld;\n",conf_max);
    }
    else
    {
        printf_s("Data \"conf_max\" from the file have not been read");
        return FileReadError;
    }

    return NoError;
}
#endif


#if __linux
int32_t Initial(void)
{
    FILE *inFile;

    inFile=fopen("Config.dat","r");//открыли файл для считывания
    if(!inFile){printf("File \"Config.dat\" no open!"); return OpenFileError;}

    //считываем данные из файла и проверяем -- получилось ли считать данные.

    //считываем температуру
    if(fscanf(inFile,"T=%lf;\n",&T))
    {
        printf("T=%lf;\n",T);
    }
    else
    {
        printf("Data \"T\" from the file have not been read");
        return FileReadError;
    }

    //считываем параметр, который отвечает за начальное состояние
    if(fscanf(inFile,"a=%lf;\n",&a))
    {
        printf("a=%lf;\n",a);
    }
    else
    {
        printf("Data \"a\" from the file have not been read");
        return FileReadError;
    }

    //считываем размер системы
    if(fscanf(inFile,"L=%ld;\n",&L))//%lld для linux заменить на ld или d?
    {
        printf("L=%ld;\n",L);
        if(L>N)
        {
            printf("L should be no more than 256");
            return RandgeError;
        }
    }
    else
    {
        printf("Data \"L\" from the file have not been read");
        return FileReadError;
    }

    //считываем максимальное количество временных шагов
    if(fscanf(inFile,"mcs_max=%ld;\n",&mcs_max))
    {
        printf("mcs_max=%ld;\n",mcs_max);
    }
    else
    {
        printf("Data \"mcs_max\" from the file have not been read");
        return FileReadError;
    }

    //считываем максимальное количество конфигураций
    if(fscanf(inFile,"conf_max=%ld;\n",&conf_max))
    {
        printf("conf_max=%ld;\n",conf_max);
    }
    else
    {
        printf("Data \"conf_max\" from the file have not been read");
        return FileReadError;
    }

    return NoError;
}
#endif

inline int8_t sosedu(int8_t inSp[N][N], size_t inI, size_t inJ)
{
    int8_t left,right,up,down;//сосед слева, справа, сверху, снизу
    if(inI==0){left=inSp[L-1][inJ];}else{left=inSp[inI-1][inJ];}
    if(inI==L-1){right=inSp[0][inJ];}else{right=inSp[inI+1][inJ];}
    if(inJ==0){down=inSp[inI][L-1];}else{down=inSp[inI][inJ-1];}
    if(inJ==L-1){up=inSp[inI][0];}else{up=inSp[inI][inJ+1];}

    return left+right+down+up;
}

inline int8_t soseduDiagonal(int8_t inSp[N][N], size_t inI, size_t inJ)
{
    int8_t left,right,up,down;//сосед слева, справа, сверху, снизу
    if(inI==0){left=L-1;}else{left=inI-1;}
    if(inI==L-1){right=0;}else{right=inI+1;}
    if(inJ==0){down=L-1;}else{down=inJ-1;}
    if(inJ==L-1){up=0;}else{up=inJ+1;}

    return inSp[left][up]+inSp[right][up]+inSp[left][down]+inSp[right][down];
}

inline double w(const double inDE, const double inT)
{
//    return exp(-inDE/T);//метрополиса
    return 1/(exp(inDE/inT)+1);//тепловая баня
}

inline double calculateMagnetic(int8_t inSp[N][N])
{
    double rez=0;
    for(size_t i=0;i<L;++i)
    {
        for(size_t j=0;j<L;++j)
        {
            rez+=inSp[i][j];
        }
    }
    return rez/(1.0*L*L);
}

inline double calculateEnergy(int8_t inSp[N][N])
{
    double rez=0;
    double J0=1;//ферромагнетик
    double J1=J0/2.0;
    for(size_t i=0;i<L;++i)
    {
        for(size_t j=0;j<L;++j)
        {
            rez+=-0.5*(J0*inSp[i][j]*sosedu(inSp,i,j)+J1*inSp[i][j]*soseduDiagonal(inSp,i,j));
        }
    }
    return rez*1.0;
}
