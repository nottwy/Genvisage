#include "rocchio.h"

void initlization(){
    FEATURE_CONSIDER_TOTAL=0; EARLY_TERMINATE = false;
    SORT_G = false; SORT_F = false;
    HORIZONTAL = true; WEIGHTED = false;
    SAMPLING = false; SAMPLING_OPT = false;
    DELIMITER = "space"; TOPK = 1000;
    DELTA = 0.05; TRASNFORM = true; PRINT = false; REST = true,HIST = false;
}

void parse_args(int argv,  char *argc[]){
  int ptr = 1;
  string cur_arg;
  
  string delim_str ="-delimiter";
  string fconsider_str = "-Fconsider";  // consider only top 1D features
  string earlyT_str = "-earlyT";  // early termination
  string sortG_str = "-sortG";  // sort genes
  string sortF_str = "-sortF";  // SORT FEATURES ACCORDING TO 1D score
  string vertical_str = "-vertical"; //horizontal traversal or vertical traversal 
  string weighted_str = "-weighted";  // weighted on pos v.s. neg
  string sampl_str = "-sampl";   // sampling
  string sampl_opt_str = "-samplOpt";  //sort the candidates according to the lower bound


  while (ptr != argv){
    cur_arg = argc[ptr];
    if (cur_arg.compare(delim_str) == 0){
        ptr++;
        DELIMITER = argc[ptr]; //-top1D 1000
    } else
    if (cur_arg.compare(fconsider_str) == 0){
        ptr++;
        cur_arg = argc[ptr];
        FEATURE_CONSIDER_TOTAL = atol(cur_arg.c_str());  //-top1D 1000
        cout << "FEATURE_CONSIDER_TOTAL =" <<FEATURE_CONSIDER_TOTAL << endl;
    } else
    if (cur_arg.compare(earlyT_str) == 0){
        EARLY_TERMINATE = true;
    } else
    if (cur_arg.compare(sortG_str) == 0){
        SORT_G = true;
    } else
    if (cur_arg.compare(sortF_str) == 0){
        SORT_F = true;
    } else
    if (cur_arg.compare(vertical_str) == 0){
        HORIZONTAL = false;
    } else
    if (cur_arg.compare(weighted_str) == 0){
        WEIGHTED = true;
    } else
    if(cur_arg.compare(sampl_str) == 0){
        SAMPLING = true;
    } else
    if(cur_arg.compare(sampl_opt_str) == 0){
        SAMPLING_OPT = true;
    }else
    if(cur_arg.compare("-matrixF") == 0){
        ptr++;
        MATRIXF = argc[ptr];
    }else
    if(cur_arg.compare("-expF") == 0){
        ptr++;
        EXPF = argc[ptr];
    }else
    if(cur_arg.compare("-expF2") == 0){
        ptr++;
        EXPF2 = argc[ptr];
    }else
    if(cur_arg.compare("-topK") == 0){
        ptr++;
        TOPK = atoi(argc[ptr]);
    }else
    if(cur_arg.compare("-delta") == 0){
        ptr++;
        DELTA = atof(argc[ptr]);
    }else
    if(cur_arg.compare("-notransform") == 0){
       TRASNFORM = false;
    }else
    if(cur_arg.compare("-print") == 0){
       PRINT = true;
       ptr++;
       printF1 = atoi(argc[ptr]);
       ptr++;
       printF2 = atoi(argc[ptr]);
       ptr++;
       buckets = atoi(argc[ptr]);
    }else
    if(cur_arg.compare("-hist") == 0){
       HIST = true;
       ptr++;
       histSize = atoi(argc[ptr]);
    }else
    if(cur_arg.compare("-pos_neg") == 0){
       REST = false;
    }
    ptr++;
  }
}

void load_feature_gene_matrix_comma(){  // matrix fille, gene is the sample
    clock_t t = clock();
    ifstream fin(MATRIXF.c_str());
    int rid =0;   // row id 

    while(fin){
        string s;
        if(!getline(fin,s)) break;
        istringstream ss(s);
        vector<string> values;
        vector<val_type > cur_f;

        int cid =0;    //column id
        while(ss){
            string value;
            if(!getline(ss,value,',')) break;   // the file's delimeter is ','
            if(rid>0 && cid!=0)      //fid = rid-1 and the first column is the feature's name
                cur_f.push_back(atof(value.c_str())); //
            else{
                if(rid ==0 && cid!=0){ //the first line describes the genes' names
                    gene_id[value.substr(1,value.size()-2)] = cid-1; // from string to id
                    //cout << cid-1 << " "<< value << " " << value.substr(1,value.size()-2)<<endl ;
                }
                else{  // rid>0 && cid =0
                    if(rid >0)//cout << rid-1 <<" " << value << endl;
                        featureName.push_back(value);  //featureName[rid-1] = value
                }
            }
            cid++;
        }
        if(rid>0){
            feature cur(cur_f,0,0);
            features.push_back(cur);
            //matrix.push_back(cur_f); //matrix[f][g]
        }
        rid++;
    }
    GENE_NUM = features[0].vals.size();  FEATURE_NUM=features.size(); MATRIX_GENE_NUM=GENE_NUM;
    //GENE_NUM = matrix[0].size();  FEATURE_NUM=matrix.size(); MATRIX_GENE_NUM=GENE_NUM;
    fprintf(stderr,"=====total number of samples: %ld==%d; total number of features: %ld==%d======\n",gene_id.size(),GENE_NUM, featureName.size(),FEATURE_NUM);
    t = clock() - t;
    fprintf (stderr, "load_feature_gene_separatedBydelimeter: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}


void load_feature_gene_matrix_space(){  // matrix fille, gene is the sample MSIGDB
    clock_t t = clock();
    ifstream fin(MATRIXF.c_str());
    //ofstream fout("test.txt");
    int rid =0, cid =0;   // row id  column id
    //int test_ccnt = 5000, test_rcnt = 50;
    string s;

    if(getline(fin,s)){
        istringstream ss(s);
        string value;
        while(ss>>value){
            gene_id[value] = cid; // from string to id
            //if(cid < test_ccnt)
            //    fout << " " << value;
            cid++;
        }
    }
    
    GENE_NUM = gene_id.size();
    cout << "GENE_NUM=" <<GENE_NUM<<";cid="<<cid<<endl;
    
    cid=0;rid++; // the second row
    vector<val_type > cur_f;
    while(fin>>s){
        if(cid ==0){
            featureName.push_back(s);
            //cout << rid << " " <<s << endl;
        }
        else
            cur_f.push_back(atof(s.c_str()));
        /*if(rid <=test_rcnt && cid < test_ccnt){
            if(cid == 0)
                fout << endl << s;
            else
                fout << " " << s;
        }*/
        cid++;
        if(cid > GENE_NUM){
            cid =0;
            rid++;
            feature cur(cur_f,0,0);
            features.push_back(cur);
            //matrix.push_back(cur_f);
            cur_f.clear();
        }
    }
        
    //cout <<"++++++++" << rid << " " << cid <<"++++++++" <<endl;
    GENE_NUM = features[0].vals.size();  FEATURE_NUM=features.size(); MATRIX_GENE_NUM=GENE_NUM;
    //GENE_NUM = matrix[0].size();  FEATURE_NUM=matrix.size(); MATRIX_GENE_NUM=GENE_NUM;
    fprintf(stderr,"----total number of samples: %ld==%d; total number of features: %ld==%d----\n",gene_id.size(),GENE_NUM, featureName.size(),FEATURE_NUM);
    t = clock() - t;
    fprintf (stderr, "load_feature_gene_matrix_space: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

/*
void load_sample_feature_matrix_space(){  // matrix fille, gene is the sample
    clock_t t = clock();
    ifstream fin(MATRIXF.c_str());
    //ofstream fout("test.txt");
    int rid =0, cid =0;   // row id  column id
    //int test_ccnt = 5000, test_rcnt = 50;
    string s;

    if(getline(fin,s)){
        istringstream ss(s);
        string value;
        while(ss>>value){
            featureName.push_back(value);
            //if(cid < test_ccnt)
            //    fout << " " << value;
            cid++;
        }
    }
    
    FEATURE_NUM = featureName.size();
    cout << "FEATURE_NUM=" <<FEATURE_NUM<<";cid="<<cid<<endl;
    matrix.resize(cid); // number of features
    
    cid=0;rid++; // the second row
    while(fin>>s){
        if(cid ==0){
            gene_id[s]=rid-1;
            //cout << rid << " " <<s << endl;
        }
        else{
            matrix[cid-1].push_back(atof(s.c_str()));
        }
        cid++;
        if(cid > FEATURE_NUM){
            cid =0;
            rid++;
        }
    }
        
    //cout <<"++++++++" << rid << " " << cid <<"++++++++" <<endl;
    GENE_NUM = matrix[0].size();  FEATURE_NUM=matrix.size(); MATRIX_GENE_NUM=GENE_NUM;
    fprintf(stderr,"----total number of samples/genes: %ld==%d; total number of features: %ld==%d----\n",gene_id.size(),GENE_NUM, featureName.size(),FEATURE_NUM);
    t = clock() - t;
    fprintf (stderr, "load_feature_gene_matrix_space: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}*/

void load_sample_feature_matrix_space(){  // matrix fille, gene is the sample LINCS
    clock_t t = clock();
    ifstream fin(MATRIXF.c_str());
    //ofstream fout("test.txt");
    int rid =0, cid =0;   // row id  column id
    //int test_ccnt = 5000, test_rcnt = 50;
    string s;

    if(getline(fin,s)){
        istringstream ss(s);
        string value;
        while(ss>>value){
            featureName.push_back(value);
            //if(cid < test_ccnt)
            //    fout << " " << value;
            cid++;
        }
    }
    
    FEATURE_NUM = featureName.size();
    cout << "FEATURE_NUM=" <<FEATURE_NUM<<";cid="<<cid<<endl;
    features.resize(cid); // number of features
    
    cid=0;rid++; // the second row
    while(fin>>s){
        if(cid ==0){
            gene_id[s]=rid-1;
            //cout << rid << " " <<s << endl;
        }
        else{
            features[cid-1].vals.push_back(atof(s.c_str()));
        }
        /*if(rid <=test_rcnt && cid < test_ccnt){
            if(cid == 0)
                fout << endl << s;
            else
                fout << " " << s;
        }*/
        cid++;
        if(cid > FEATURE_NUM){
            cid =0;
            rid++;
        }
    }
        
    //cout <<"++++++++" << rid << " " << cid <<"++++++++" <<endl;
    GENE_NUM = features[0].vals.size();  FEATURE_NUM=features.size(); MATRIX_GENE_NUM=GENE_NUM;
    fprintf(stderr,"----total number of samples/genes: %ld==%d; total number of features: %ld==%d----\n",gene_id.size(),GENE_NUM, featureName.size(),FEATURE_NUM);
    t = clock() - t;
    fprintf (stderr, "load_feature_gene_matrix_space: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}


void histogram_feature_sample_matrix_comma(){  // matrix fille, gene is the sample 
    clock_t t = clock();
    ifstream fin(MATRIXF.c_str());
    int rid =0;   // row id 
    std::vector<val_type> vals;
    
    while(fin){
        string s;
        if(!getline(fin,s)) break;
        istringstream ss(s);
        vector<string> values;
        vector<val_type > cur_f;

        int cid =0;    //column id
        while(ss){
            string value;
            if(!getline(ss,value,',')) break;   // the file's delimeter is ','
            if(rid>0 && cid!=0)      //fid = rid-1 and the first column is the feature's name
                vals.push_back(atof(value.c_str()));
            else{
                if(rid ==0 && cid!=0){ //the first line describes the genes' names
                    gene_id[value.substr(1,value.size()-2)] = cid-1; // from string to id
                    //cout << cid-1 << " "<< value << " " << value.substr(1,value.size()-2)<<endl ;
                }
                else{  // rid>0 && cid =0
                    if(rid >0)//cout << rid-1 <<" " << value << endl;
                        featureName.push_back(value);  //featureName[rid-1] = value
                }
            }
            cid++;
        }
        rid++;
    }
    sort(vals.begin(),vals.end());
    
    FILE* fhist = fopen("histgram_freq.txt","w");
    int cur = 0,step = vals.size()/histSize;
    cout << "frequent: step size " << step <<endl;
    for(int i=0; i<histSize; i++){
        if(i == histSize-1){
            fprintf(fhist, "%f \t %f \t\t %d \n", vals[cur], vals[vals.size()-1], vals.size()-cur);
            break;
        }
        fprintf(fhist, "%f \t %f \t\t %d \n", vals[cur], vals[cur+step], step);
        cur += step;
    }
    fclose(fhist);
    
    FILE* whist = fopen("histgram_width.txt","w");
    double val_step = (vals[vals.size()-1] - vals[0])/histSize;
    cout << "width: step size " << val_step <<endl;
    double cur_boundary = vals[0];
    int cnt = 0, j=0;
    
    for(int i=0; i<histSize; i++){
        
        while(j < vals.size() && (vals[j] < cur_boundary + val_step || j == vals.size()-1))
        {
            cnt++; j++;
        }
        fprintf(whist, "%f \t %f \t\t %d \n", cur_boundary, cur_boundary + val_step, cnt);
        cur_boundary += val_step;
        cnt = 1;
    }
    fclose(fhist);
    
    fprintf (stderr, "Histogram: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

void histogram_feature_sample_matrix_space(){  // matrix fille, gene is the sample 
    clock_t t = clock();
    ifstream fin(MATRIXF.c_str());
    int rid =0, cid =0;   // row id  column id
    string s;

    if(getline(fin,s)){
        istringstream ss(s);
        string value;
        while(ss>>value){
            gene_id[value] = cid; // from string to id
            cid++;
        }
    }
    
    GENE_NUM = gene_id.size();
    cout << "GENE_NUM=" <<GENE_NUM<<";cid="<<cid<<endl;
    
    std::vector<val_type> vals;
    cid=0;rid++; // the second row
    while(fin>>s){
        if(cid ==0){
            featureName.push_back(s);
        }
        else
            vals.push_back(atof(s.c_str()));
        cid++;
        if(cid > GENE_NUM){
            cid =0;
            rid++;
        }
    }
    sort(vals.begin(),vals.end());
    
    FILE* fhist = fopen("histgram_freq.txt","w");
    int cur = 0,step = vals.size()/histSize;
    cout << "frequent: step size " << step <<endl;
    for(int i=0; i<histSize; i++){
        if(i == histSize-1){
            fprintf(fhist, "%f \t %f \t\t %d \n", vals[cur], vals[vals.size()-1], vals.size()-cur);
            break;
        }
        fprintf(fhist, "%f \t %f \t\t %d \n", vals[cur], vals[cur+step], step);
        cur += step;
    }
    fclose(fhist);
    
    FILE* whist = fopen("histgram_width.txt","w");
    double val_step = (vals[vals.size()-1] - vals[0])/histSize;
    cout << "width: step size " << val_step <<endl;
    double cur_boundary = vals[0];
    int cnt = 0, j=0;
    
    for(int i=0; i<histSize; i++){
        
        while(j < vals.size() && (vals[j] < cur_boundary + val_step || j == vals.size()-1))
        {
            cnt++; j++;
        }
        fprintf(whist, "%f \t %f \t\t %d \n", cur_boundary, cur_boundary + val_step, cnt);
        cur_boundary += val_step;
        cnt = 1;
    }
    fclose(fhist);
    
    fprintf (stderr, "Histogram: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}


void histogram_sample_feature_matrix_space(){  // matrix fille, gene is the sample
    clock_t t = clock();
    ifstream fin(MATRIXF.c_str());
    int rid =0, cid =0;   // row id  column id
    string s;

    if(getline(fin,s)){
        istringstream ss(s);
        string value;
        while(ss>>value){
            cid++;
        }
    }
    
    FEATURE_NUM = cid;
    cout << "FEATURE_NUM=" <<FEATURE_NUM<<";cid="<<cid<<endl;
    
    std::vector<val_type> vals;
    cid=0;rid++; // the second row
    while(fin>>s){
        if(cid ==0){
            gene_id[s]=rid-1;
        }
        else{
            vals.push_back(atof(s.c_str()));
        }
        cid++;
        if(cid > FEATURE_NUM){
            cid =0;
            rid++;
        }
    }
    sort(vals.begin(),vals.end());
    
    FILE* fhist = fopen("histgram_freq.txt","w");
    int cur = 0,step = vals.size()/histSize;
    cout << "frequent: step size " << step <<endl;
    for(int i=0; i<histSize; i++){
        if(i == histSize-1){
            fprintf(fhist, "%f \t %f \t\t %d \n", vals[cur], vals[vals.size()-1], vals.size()-cur);
            break;
        }
        fprintf(fhist, "%f \t %f \t\t %d \n", vals[cur], vals[cur+step], step);
        cur += step;
    }
    fclose(fhist);
    
    FILE* whist = fopen("histgram_width.txt","w");
    double val_step = (vals[vals.size()-1] - vals[0])/histSize;
    cout << "width: step size " << val_step <<endl;
    double cur_boundary = vals[0];
    int cnt = 0;
    for(int i=0; i<vals.size();i++){
        if(i == vals.size()-1){
            fprintf(whist, "%f \t %f \t\t %d \n", cur_boundary, vals[vals.size()-1], cnt);
            break;
        }
        if(vals[i] < cur_boundary + val_step){
            cnt++;
        }else{
            fprintf(whist, "%f \t %f \t\t %d \n", cur_boundary, cur_boundary + val_step, cnt);
            cur_boundary += val_step;
            cnt = 1;
        }
    }
    fclose(fhist);
    
    fprintf (stderr, "Histogram: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}


/*void load_feature_gene_matrix_space(){  // matrix fille, gene is the sample
    clock_t t = clock();
    FEATURE_NUM = 21202; GENE_NUM=22209; 
    FILE * fin = fopen(MATRIXF.c_str(),"r");
    matrix.resize(FEATURE_NUM);

    char name[200];
    for(int gid=0;gid<GENE_NUM;gid++) {
        int x=fscanf(fin,"%s", name);  // gene name 
        cout << gid << " " << name << endl;
        string gene_name(name);
        gene_id[gene_name] = gid;  // gene id
    }

    for(int f=0; f<FEATURE_NUM; f++){
        int x = fscanf(fin,"%s ", name);
        featureName.push_back(name);
        cout << f<< " " << name  << endl;
        val_type entry;
        for(int g=0; g<GENE_NUM;g++){
            x=fscanf(fin,"%lf ", &entry);
            //if(f==0)
            //fprintf(stderr,"=====%d,%lf======\n",g,entry);
            matrix[f].push_back(entry);
        }
    }
    
    fprintf(stderr,"----total number of samples: %ld==%d; total number of features: %ld==%d----\n",gene_id.size(),GENE_NUM, featureName.size(),FEATURE_NUM);
    t = clock() - t;
    fprintf (stderr, "load_feature_gene_matrix_space: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}*/
void load_exp_pos_neg(){  //not all samples are selected
    clock_t t = clock();
    FILE *fgene1 = fopen(EXPF.c_str(),"r");
    FILE *fgene2 = fopen(EXPF2.c_str(),"r");

    char name[200];
    pos_num=0;
    genes.resize(GENE_NUM, 0); //positive or negative
    if(fgene1 != NULL){
        while(fscanf(fgene1,"%s ", name)!=EOF){
            string gene_name(name);
            map<string, int>::iterator it = gene_id.find(gene_name);
            if(it != gene_id.end()){
                genes[it->second]=1;
                pos_num++;
            }
        }
    }
    fprintf(stderr,"=====number of positive genes:%d;======\n",pos_num);
    int neg_num =0;
    if(fgene2 != NULL){
        while(fscanf(fgene2,"%s ", name)!=EOF){
            string gene_name(name);
            map<string, int>::iterator it = gene_id.find(gene_name);
            if(it != gene_id.end()){
                genes[it->second]=-1;
                neg_num++;
            }
        }
    }
    fprintf(stderr,"=====number of negative genes:%d;======\n",neg_num);
    
    features.resize(FEATURE_NUM);
    for(int f=0;f<FEATURE_NUM;f++){
        features[f].fid = f;
        vector<val_type > cur = features[f].vals;
        int gg = 0;
        for(int g=0; g<GENE_NUM;g++){
            if(genes[g]==1){//pos genes
                features[f].vals[gg] = cur[g];
                gg++;
                //features[f].vals.push_back(matrix[f][g]);
            }
        }
        for(int g=0; g<GENE_NUM;g++){
            if(genes[g]==-1){//neg genes
                features[f].vals[gg] = cur[g];
                gg++;
                //features[f].vals.push_back(matrix[f][g]);
            }
        }
    }
    GENE_NUM = pos_num+ neg_num;
    fprintf(stderr,"------total number of genes:%ld==%d-------\n",features[0].vals.size(), GENE_NUM);
    fprintf(stderr,"------done with extracting submatrix for this experiment-------\n");
    t = clock() - t;
    fprintf (stderr, "load exp: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}


/*
void load_exp_rest(){
    clock_t t = clock();
    FILE *fgene1 = fopen(EXPF.c_str(),"r");
    FILE *fgene2 = fopen(EXPF2.c_str(),"r");

    char name[200];
    pos_num=0;
    genes.resize(GENE_NUM, -1); //positive or negative
    if(fgene1 != NULL){
        while(fscanf(fgene1,"%s ", name)!=EOF){
            string gene_name(name);
            map<string, int>::iterator it = gene_id.find(gene_name);
            if(it != gene_id.end()){
                genes[it->second]=1;
                pos_num++;
            }
        }
    }
    fprintf(stderr,"=====number of positive genes:%d;======\n",pos_num);
    if(fgene2 != NULL){
        while(fscanf(fgene2,"%s ", name)!=EOF){
            string gene_name(name);
            map<string, int>::iterator it = gene_id.find(gene_name);
            if(it != gene_id.end()){
                genes[it->second]=1;
                pos_num++;
            }
        }
    }
    fprintf(stderr,"=====number of positive genes:%d;======\n",pos_num);
    fprintf(stderr,"=====FEATURE_NUM:%d;GENE_NUM:%d======\n",FEATURE_NUM,GENE_NUM);
    
    features.resize(FEATURE_NUM);
    for(int f=0;f<FEATURE_NUM;f++){
        features[f].fid = f;
        for(int g=0; g<GENE_NUM;g++){
            if(genes[g]==1){//pos genes
                features[f].vals.push_back(matrix[f][g]);
            }
        }
        for(int g=0; g<GENE_NUM;g++){
            if(genes[g]==-1){//neg genes
                features[f].vals.push_back(matrix[f][g]);
            }
        }
    }
    
    fprintf(stderr,"------total number of genes:%ld==%d-------\n",features[0].vals.size(), GENE_NUM);
    fprintf(stderr,"------done with extracting submatrix for this experiment-------\n");
    t = clock() - t;
    fprintf (stderr, "load exp: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}*/


/// after load_exp_rest(), "features[f].vals" are ordered by positive samples and then negative samples;
void load_exp_rest(){
    clock_t t = clock();
    FILE *fgene1 = fopen(EXPF.c_str(),"r");
    FILE *fgene2 = fopen(EXPF2.c_str(),"r");

    char name[200];
    pos_num=0;
    genes.resize(GENE_NUM, -1); //positive or negative
    if(fgene1 != NULL){
        while(fscanf(fgene1,"%s ", name)!=EOF){
            string gene_name(name);
            map<string, int>::iterator it = gene_id.find(gene_name);
            if(it != gene_id.end()){
                genes[it->second]=1;
                pos_num++;
            }
        }
    }
    fprintf(stderr,"=====number of positive genes:%d;======\n",pos_num);
    if(fgene2 != NULL){
        while(fscanf(fgene2,"%s ", name)!=EOF){
            string gene_name(name);
            map<string, int>::iterator it = gene_id.find(gene_name);
            if(it != gene_id.end() && genes[it->second]!=1){
                genes[it->second]=1;
                pos_num++;
            }
        }
    }
    fprintf(stderr,"=====number of positive genes:%d;======\n",pos_num);
    fprintf(stderr,"=====FEATURE_NUM:%d;GENE_NUM:%d======\n",FEATURE_NUM,GENE_NUM);
    
    features.resize(FEATURE_NUM);
    for(int f=0;f<FEATURE_NUM;f++){
        features[f].fid = f;
        vector<val_type > cur = features[f].vals;
        int gg = 0;
        for(int g=0; g<GENE_NUM;g++){
            if(genes[g]==1){//pos genes
                features[f].vals[gg] = cur[g];
                gg++;
            }
        }
        for(int g=0; g<GENE_NUM;g++){
            if(genes[g]==-1){//neg genes
                features[f].vals[gg] = cur[g];
                gg++;
            }
        }
    }
    
    fprintf(stderr,"------total number of genes:%ld==%d-------\n",features[0].vals.size(), GENE_NUM);
    fprintf(stderr,"------done with extracting submatrix for this experiment-------\n");
    t = clock() - t;
    fprintf (stderr, "load exp: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}
     
/*transformation and get top features*/
void transformation(){
    clock_t t = clock();
    ff.resize(FEATURE_NUM);
    cout << "pos_num="<<pos_num<<";GENE_NUM="<<GENE_NUM<<endl;
    priority_queue<top_1d> top_1d_features;
    for(int f=0; f<FEATURE_NUM; f++){
        if(features[f].vals.size()!=GENE_NUM)
            cout << "*******ERROR******"<<f<<" size:"<< features[f].vals.size() << endl;
        vector<val_type > pos_genes(features[f].vals.begin(), features[f].vals.begin()+pos_num);
        vector<val_type > neg_genes(features[f].vals.begin()+pos_num, features[f].vals.end());

        sort(pos_genes.begin(), pos_genes.end());
        val_type median_p =pos_genes[pos_genes.size()/2+1];  //x+
        sort(neg_genes.begin(), neg_genes.end());
        val_type median_n =neg_genes[neg_genes.size()/2+1]; //x-     
        val_type w = median_p - median_n;  //(x+ - x-)
        val_type intercept = -(median_p*median_p-median_n*median_n)/2;  //-(x+^2- x-^2)/2
        
        int corP=0,corN=0; //correct positive genes & correct negative genes
        for(int i=0; i<pos_num;i++){
            val_type tmp=features[f].vals[i] * w + intercept; //CHANGE: int TO val_type
            features[f].vals[i] = tmp; 

            if(features[f].vals[i]>0){
                corP++;
            }
        }
        for(int i=pos_num; i<GENE_NUM; i++){
            val_type tmp= -(features[f].vals[i] * w + intercept); ////CHANGE: int TO val_type
            features[f].vals[i] = tmp;
            if(features[f].vals[i]>0){
                corN++;
            }
        }
   
        features[f].correct = corP*pos_weight+corN;  //weighted!!!
        top_1d cur_f(f,corP*pos_weight+corN,corP,corN);  // can be omitted 
        top_1d_features.push(cur_f);  // can be omitted
        ff[f] = features[f];
    }
    fprintf(stderr,"------done with transformation-------\n");
    t = clock() - t;
    fprintf (stderr, "transformation: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    
    FILE *ftop1d = fopen("FTOP1d.txt","a");
    //correfprintf(ftop1d,"fid\tcorrect\terror\tpos_correct\tneg_correct\n");
    for(int i=0; i<1000; i++){
        top_1d cur = top_1d_features.top();
        fprintf(ftop1d,"%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",cur.fid, featureName[cur.fid].c_str(), cur.correct,pos_num*pos_weight+(GENE_NUM-pos_num)-cur.correct,cur.correct_p,cur.correct_n,pos_num-cur.correct_p,GENE_NUM-pos_num-cur.correct_n);
        if(i<10)
            fprintf(stderr,"%d correct in top-%d feature %d \n",cur.correct,i,cur.fid);
        top_1d_features.pop();
    }
    int rank = 1000;
    /*while(!top_1d_features.empty()){
        top_1d cur = top_1d_features.top();
        top_1d_features.pop();
        if(cur.fid==15398)
            fprintf(stderr,"%d correct in top-%d feature %d \n",cur.correct,rank,cur.fid);
        rank++;
    }*/
    fclose(ftop1d);
}


void no_transformation_baseline(FILE* ftime){
    clock_t t = clock();
    
    top_2d best[TOPK];
    for(int i=0;i<TOPK; i++)
        best[i].miss=2*GENE_NUM;
    /*each feature pair*/

    long total_checked_g=0;
    for(int f1=0;f1<FEATURE_NUM;f1++){
        vector<val_type > pos_genes_f1(features[f1].vals.begin(), features[f1].vals.begin()+pos_num);
        vector<val_type > neg_genes_f1(features[f1].vals.begin()+pos_num, features[f1].vals.end());
        //median of positive 
        sort(pos_genes_f1.begin(), pos_genes_f1.end());
        val_type median_p_f1 =pos_genes_f1[pos_genes_f1.size()/2+1];  //x+
        sort(neg_genes_f1.begin(), neg_genes_f1.end());
        val_type median_n_f1 =neg_genes_f1[neg_genes_f1.size()/2+1]; //x-
        val_type w_f1 = median_p_f1 - median_n_f1;  //(x+ - x-)
        for(int f2=f1+1; f2<FEATURE_NUM;f2++){
            vector<val_type > pos_genes_f2(features[f2].vals.begin(), features[f2].vals.begin()+pos_num);
            vector<val_type > neg_genes_f2(features[f2].vals.begin()+pos_num, features[f2].vals.end());
        //median of positive 
            sort(pos_genes_f2.begin(), pos_genes_f2.end());
            val_type median_p_f2 =pos_genes_f2[pos_genes_f2.size()/2+1];  //x+
            sort(neg_genes_f2.begin(), neg_genes_f2.end());
            val_type median_n_f2 =neg_genes_f2[neg_genes_f2.size()/2+1]; //x-
        
            val_type w_f2 = median_p_f2 - median_n_f2;  //(y+ - y-)
            val_type intercept = -(median_p_f1*median_p_f1-median_n_f1*median_n_f1)/2-(median_p_f2*median_p_f2-median_n_f2*median_n_f2)/2;  //-(x+^2- x-^2)/2
            //cout << "w_f1 =" << w_f1 <<"; w_f2="<<w_f2 <<";intercept="<<intercept<< endl;
            int g=0, missed_p=0, missed_n=0,missed =0;
            for(; g<pos_num; g++){
                val_type result = w_f1*features[f1].vals[g]+w_f2*features[f2].vals[g]+intercept;
                if(result <= 0){ // label is positive
                    missed_p++;
                }
            }
            for(;g<GENE_NUM;g++){
                val_type result = w_f1*features[f1].vals[g]+w_f2*features[f2].vals[g]+intercept;
                if(result >= 0){ //!!! label is negative
                    missed_n++;
                }
            }
            //cout << "missed_p = " << missed_p << ";missed_n = " << missed_n << endl; 
            missed=pos_weight*missed_p+missed_n;
            
            if(missed < best[0].miss){
                top_2d curFpair(f1,f2,missed,missed_p,missed_n);
                best[0] = curFpair;
                sort(best,best+TOPK);
            }
            total_checked_g+=g;
        }
        if(f1%1000==0)
            cout << "f1=" << f1 << "; total_checked_g=" <<total_checked_g<<endl;
    }

    t = clock() - t;
    int per_checked_g = 2*total_checked_g/(FEATURE_NUM-1)/FEATURE_NUM;
    cout << "total_checked_g="<<total_checked_g<<"; per_checked_g="<< per_checked_g << endl;
    
    fprintf (stderr, "******best feature pair (%d,%d:%d) (horizontal)*****\n",features[best[TOPK-1].f1].fid,features[best[TOPK-1].f2].fid,best[TOPK-1].miss);
    fprintf (stderr, "no_transformation_baseline find best feature pair: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    fprintf(ftime,"%ld\t%d\t%f\n",total_checked_g,per_checked_g,((float)t)/CLOCKS_PER_SEC);
    
    FILE *no_transform_file = fopen("FP_no_transform.txt","a");
    for(int i=TOPK-1;i>=0;i--){
        int f1_tmp= features[best[i].f1].fid, f2_tmp=features[best[i].f2].fid;
        int f1 = min(f1_tmp,f2_tmp), f2 = max(f1_tmp,f2_tmp);
        fprintf(no_transform_file,"%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", f1,f2,featureName[f1].c_str(),featureName[f2].c_str(),pos_num*pos_weight+(GENE_NUM-pos_num)-best[i].miss, best[i].miss,pos_num-best[i].miss_p, GENE_NUM-pos_num-best[i].miss_n, best[i].miss_p,best[i].miss_n);
    }
    fclose(no_transform_file);
}



void sort_genes_weighted(){
    clock_t t = clock();
    ff.clear();
    ff.resize(FEATURE_NUM);
    vector<topG > gene_correct_p,gene_correct_n;  // use top_1d (id feature)-- fid->gid
    for(int g=0;g<GENE_NUM;g++){
        int g_corrects = 0;
        for(int f=0;f<FEATURE_NUM;f++){
            if(features[f].vals[g]>0)
                g_corrects++;
        }
        topG tmp(g,g_corrects);
        if(g<pos_num)
            gene_correct_p.push_back(tmp);
        else
            gene_correct_n.push_back(tmp);
    }
    sort(gene_correct_p.begin(),gene_correct_p.end());
    cout << gene_correct_p[0].correct << "---test--" <<gene_correct_p[gene_correct_p.size()/2].correct<< "---test--" <<gene_correct_p[gene_correct_p.size()-1].correct<< endl;
    sort(gene_correct_n.begin(),gene_correct_n.end());
    cout << gene_correct_n[0].correct << "---test--" <<gene_correct_n[gene_correct_n.size()/2].correct<< "---test--" <<gene_correct_n[gene_correct_n.size()-1].correct<< endl;
    FILE* testf = fopen("gene_reordering.txt", "w");

    cout << "pos_num=" << pos_num << ";GENE_NUM=" <<GENE_NUM << endl;
    for(int f=0; f<FEATURE_NUM;f++){
        ff[f].correct = features[f].correct;
        ff[f].fid = features[f].fid;
        for(int g=0; g<pos_num;g++){
            ff[f].vals.push_back(features[f].vals[gene_correct_p[g].gid]);  //gene_correct[g].fid: gid of gth gene
            if(f==0)
                fprintf(testf,"%d \t %d \t %d\n",g, gene_correct_p[g].gid, gene_correct_p[g].correct);

        }
        for(int g=0; g<GENE_NUM-pos_num;g++){
            ff[f].vals.push_back(features[f].vals[gene_correct_n[g].gid]);  //gene_correct[g].fid: gid of gth gene
            if(f==0)
                fprintf(testf,"%d \t %d \t %d\n",g, gene_correct_n[g].gid,  gene_correct_n[g].correct);
        }
        if(ff[f].vals.size() > GENE_NUM)
            fprintf(stderr,"%d \t %ld \t %d\n",f, ff[f].vals.size(),  GENE_NUM);
    }
    fclose(testf);
    t = clock() - t;
    fprintf (stderr, "sort_genes_weighted : %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

void sort_genes(){
    clock_t t = clock();
    ff.clear();
    ff.resize(FEATURE_NUM);
    vector<topG > gene_correct;  // use top_1d (id feature)-- fid->gid
    for(int g=0;g<GENE_NUM;g++){
        int g_corrects = 0;
        for(int f=0;f<FEATURE_NUM;f++){
            if(features[f].vals[g]>0)
                g_corrects++;
        }
        topG tmp(g,g_corrects);
        gene_correct.push_back(tmp);
    }
    sort(gene_correct.begin(),gene_correct.end());
    
    for(int f=0; f<FEATURE_NUM;f++){
        ff[f].correct = features[f].correct;
        ff[f].fid = features[f].fid;
        for(int g=0; g<GENE_NUM;g++){
            ff[f].vals.push_back(features[f].vals[gene_correct[g].gid]);  //gene_correct[g].fid: gid of gth gene
        }
        if(ff[f].vals.size() > GENE_NUM)
            fprintf(stderr,"***********%d \t %ld \t %d**********\n",f, ff[f].vals.size(),  GENE_NUM);
    }
    t = clock() - t;
    fprintf (stderr, "sort_genes : %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}

void baseline(FILE *ftime){
    clock_t t = clock();
    top_2d best[TOPK];
    for(int i=0;i<TOPK; i++)
        best[i].miss=2*GENE_NUM;

    long total_checked_g=0;
    for(int f1=0;f1<FEATURE_CONSIDER;f1++){
        for(int f2=f1+1; f2<FEATURE_NUM;f2++){
            int g=0, missed_p=0, missed_n=0,missed =0;
            for(; g<pos_num; g++){
                if(ff[f1].vals[g]+ff[f2].vals[g] <= 0){
                    missed_p++;
                }
            }
            for(;g<GENE_NUM;g++){
                if(ff[f1].vals[g]+ff[f2].vals[g] <= 0){
                    missed_n++;
                }
            }
            missed=pos_weight*missed_p+missed_n;
                    
            if(missed < best[0].miss){
                top_2d curFpair(f1,f2,missed,missed_p,missed_n);
                best[0] = curFpair;
                sort(best,best+TOPK);
            }
            total_checked_g+=g;
        }
    }

    t = clock() - t;
    int per_checked_g = 2*total_checked_g/(2*FEATURE_NUM-FEATURE_CONSIDER-1)/(FEATURE_CONSIDER);
    cout << "total_checked_g="<<total_checked_g<<"; per_checked_g="<< per_checked_g << endl;
    
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered(horizontal)*****\n",ff[best[TOPK-1].f1].fid,ff[best[TOPK-1].f2].fid,best[TOPK-1].miss, FEATURE_CONSIDER);
    fprintf (stderr, "baseline find best feature pair: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    fprintf(ftime,"%ld\t%d\t%f\n",total_checked_g,per_checked_g,((float)t)/CLOCKS_PER_SEC);
    
    
    FILE *rocchio_file = fopen("FP_baseline.txt","a");
    for(int i=TOPK-1;i>=0;i--){
        int f1_tmp= ff[best[i].f1].fid, f2_tmp=ff[best[i].f2].fid;
        int f1 = min(f1_tmp,f2_tmp), f2 = max(f1_tmp,f2_tmp);
        fprintf(rocchio_file,"%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", f1,f2,featureName[f1].c_str(),featureName[f2].c_str(),pos_num*pos_weight+(GENE_NUM-pos_num)-best[i].miss, best[i].miss,pos_num-best[i].miss_p, GENE_NUM-pos_num-best[i].miss_n, best[i].miss_p,best[i].miss_n);
    }
    fclose(rocchio_file);
}


void baseline_vertical(FILE *ftime){
    clock_t t = clock();
    top_2d best[TOPK];
    for(int i=0;i<TOPK; i++)
        best[i].miss=2*GENE_NUM;

    long total_checked_g=0;
    for(int f1=1;f1<FEATURE_CONSIDER;f1++){
        for(int f2=0; f2<f1;f2++){
            int g=0, missed_p=0, missed_n=0,missed =0;
            for(; g<pos_num; g++){
                if(ff[f1].vals[g]+ff[f2].vals[g] <= 0){
                    missed_p++;
                }
            }
            for(;g<GENE_NUM;g++){
                if(ff[f1].vals[g]+ff[f2].vals[g] <= 0){
                    missed_n++;
                }
            }
            missed=pos_weight*missed_p+missed_n;
                    
            if(missed < best[0].miss){
                top_2d curFpair(f1,f2,missed,missed_p,missed_n);
                best[0] = curFpair;
                sort(best,best+TOPK);
            }
            total_checked_g+=g;
        }
    }

    t = clock() - t;
    int per_checked_g = 2*total_checked_g/(FEATURE_CONSIDER-1)/(FEATURE_CONSIDER);
    cout << "total_checked_g="<<total_checked_g<<"; per_checked_g="<< per_checked_g << endl;
    
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered(horizontal)*****\n",ff[best[TOPK-1].f1].fid,ff[best[TOPK-1].f2].fid,best[TOPK-1].miss, FEATURE_CONSIDER);
    fprintf (stderr, "baseline_vertical find best feature pair: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    fprintf(ftime,"%ld\t%d\t%f\n",total_checked_g,per_checked_g,((float)t)/CLOCKS_PER_SEC);
    
    
    FILE *rocchio_file = fopen("FP_baseline.txt","a");
    for(int i=TOPK-1;i>=0;i--){
        int f1_tmp= ff[best[i].f1].fid, f2_tmp=ff[best[i].f2].fid;
        int f1 = min(f1_tmp,f2_tmp), f2 = max(f1_tmp,f2_tmp);
        fprintf(rocchio_file,"%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", f1,f2,featureName[f1].c_str(),featureName[f2].c_str(),pos_num*pos_weight+(GENE_NUM-pos_num)-best[i].miss, best[i].miss,pos_num-best[i].miss_p, GENE_NUM-pos_num-best[i].miss_n, best[i].miss_p,best[i].miss_n);
    }
    fclose(rocchio_file);
}



void early_term_horizontal(FILE *ftime){
    clock_t t = clock();
    top_2d best[TOPK];
    for(int i=0;i<TOPK; i++)
        best[i].miss=2*GENE_NUM;

    long total_checked_g=0;
    for(int f1=0;f1<FEATURE_CONSIDER;f1++){
        for(int f2=f1+1; f2<FEATURE_NUM;f2++){
            int g=0, missed_p=0, missed_n=0,missed =0;
            for(; g<pos_num; g++){
                if(ff[f1].vals[g]+ff[f2].vals[g] <= 0){
                    missed_p++;
                    missed+=pos_weight;
                }
                if(missed >= best[0].miss)  // need or not??
                    break;
            }
            for(;g<GENE_NUM;g++){
                if(ff[f1].vals[g]+ff[f2].vals[g] <= 0){
                    missed_n++;
                    missed++;
                }
                if(missed >= best[0].miss)
                    break;
            }
                    
            if(missed < best[0].miss){
                top_2d curFpair(f1,f2,missed,missed_p,missed_n);
                best[0] = curFpair;
                sort(best,best+TOPK);
            }
            total_checked_g+=g;
        }
    }

    t = clock() - t;
    int per_checked_g = 2*total_checked_g/(2*FEATURE_NUM-FEATURE_CONSIDER-1)/(FEATURE_CONSIDER);
    cout << "total_checked_g="<<total_checked_g<<"; per_checked_g="<< per_checked_g << endl;
    
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered(horizontal)*****\n",ff[best[TOPK-1].f1].fid,ff[best[TOPK-1].f2].fid,best[TOPK-1].miss, FEATURE_CONSIDER);
    fprintf (stderr, "early_term_horizontal find best feature pair: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    fprintf(ftime,"%ld\t%d\t%f\n",total_checked_g,per_checked_g,((float)t)/CLOCKS_PER_SEC);
    
    
    FILE *rocchio_file = fopen("FP_nosampling.txt","a");
    for(int i=TOPK-1;i>=0;i--){
        int f1_tmp= ff[best[i].f1].fid, f2_tmp=ff[best[i].f2].fid;
        int f1 = min(f1_tmp,f2_tmp), f2 = max(f1_tmp,f2_tmp);
        fprintf(rocchio_file,"%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", f1,f2,featureName[f1].c_str(),featureName[f2].c_str(),pos_num*pos_weight+(GENE_NUM-pos_num)-best[i].miss, best[i].miss,pos_num-best[i].miss_p, GENE_NUM-pos_num-best[i].miss_n, best[i].miss_p,best[i].miss_n);
    }
    fclose(rocchio_file);
}


void early_term_vertical(FILE *ftime){
    clock_t t = clock();
    top_2d best[TOPK];
    for(int i=0;i<TOPK; i++)
        best[i].miss=2*GENE_NUM;

    long total_checked_g=0;
    for(int f1=1;f1<FEATURE_CONSIDER;f1++){
        for(int f2=0; f2<f1;f2++){
            int g=0, missed_p=0, missed_n=0,missed =0;
            for(; g<pos_num; g++){
                if(ff[f1].vals[g]+ff[f2].vals[g] <= 0){
                    missed_p++;
                    missed+=pos_weight;
                }
                if(missed >= best[0].miss)  // need or not??
                    break;
            }
            for(;g<GENE_NUM;g++){
                if(ff[f1].vals[g]+ff[f2].vals[g] <= 0){
                    missed_n++;
                    missed++;
                }
                if(missed >= best[0].miss)
                    break;
            }
                    
            if(missed < best[0].miss){
                top_2d curFpair(f1,f2,missed,missed_p,missed_n);
                best[0] = curFpair;
                sort(best,best+TOPK);
            }
            total_checked_g+=g;
        }
    }

    t = clock() - t;
    int per_checked_g = 2*total_checked_g/(FEATURE_CONSIDER-1)/(FEATURE_CONSIDER);
    cout << "total_checked_g="<<total_checked_g<<"; per_checked_g="<< per_checked_g << endl;
    
    fprintf (stderr, "******best feature pair (%d,%d:%d) among %d features considered(vertical)*****\n",ff[best[TOPK-1].f1].fid,ff[best[TOPK-1].f2].fid,best[TOPK-1].miss, FEATURE_CONSIDER);
    fprintf (stderr, "early_term_vertical find best feature pair: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    fprintf(ftime,"%ld\t%d\t%f\n",total_checked_g,per_checked_g,((float)t)/CLOCKS_PER_SEC);
    
    
    FILE *rocchio_file = fopen("FP_nosampling.txt","a");
    for(int i=TOPK-1;i>=0;i--){
        int f1_tmp= ff[best[i].f1].fid, f2_tmp=ff[best[i].f2].fid;
        int f1 = min(f1_tmp,f2_tmp), f2 = max(f1_tmp,f2_tmp);
        fprintf(rocchio_file,"%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", f1,f2,featureName[f1].c_str(),featureName[f2].c_str(),pos_num*pos_weight+(GENE_NUM-pos_num)-best[i].miss, best[i].miss,pos_num-best[i].miss_p, GENE_NUM-pos_num-best[i].miss_n, best[i].miss_p,best[i].miss_n);

    }
    fclose(rocchio_file);
}

void prepare_samples(){
    clock_t t = clock();
    int sampleCnt = 1/(DELTA*DELTA);
    cout << "sampleCnt=" << sampleCnt <<endl;
    sample_f.resize(FEATURE_NUM);

    if(WEIGHTED){
        for(int f=0; f< FEATURE_NUM; f++){
            for(int g=0; g<pos_num; g++){
                sample_f[f].push_back(ff[f].vals[g]);    
            }
            for(int g=0; g<sampleCnt;g++){
                int randVal = rand()%(GENE_NUM-pos_num)+pos_num; 
                sample_f[f].push_back(ff[f].vals[randVal]);
            }
        }        
    }
    else{
        for(int f=0; f< FEATURE_NUM; f++){
            for(int g=0; g<sampleCnt; g++){
                int randVal = rand()%GENE_NUM; 
                sample_f[f].push_back(ff[f].vals[randVal]);
            }
        }   
    }
    t = clock() - t;
    fprintf (stderr, "prepare_sample_data: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
}


vector<top_2d > sampling_candidate_generation_horizontal(FILE* ftime){
    clock_t t = clock();

    top_2d best[TOPK];
    vector<top_2d > candidates;
    for(int i=0;i<TOPK; i++){
        best[i].miss=2*GENE_NUM;
    }
    long total_checked_g=0,FP_cnt=0;
    int sample_pCnt = (WEIGHTED==true)? pos_num:0; 
    int sample_nCnt = 1/(DELTA*DELTA);
    int sample_allCnt = sample_nCnt+sample_pCnt; 
    int interval = (WEIGHTED==true)?(GENE_NUM-pos_num)*DELTA:GENE_NUM*DELTA;  // confidence interval
    int sample_w = (WEIGHTED == true)? (GENE_NUM-pos_num)/sample_nCnt: GENE_NUM/sample_nCnt;

    //WEIGHTED: (GENE_NUM-pos_num)*DELTA; UNWEIGHTED:GENE_NUM*DELTA
    cout<<"sample_pCnt="<<sample_pCnt <<";sample_allCnt="<<sample_allCnt<<endl;
    cout << "confidence interval = " << interval << endl;
    int total_fp =0;
    for(int f1=0;f1<FEATURE_CONSIDER;f1++){
        for(int f2=f1+1; f2<FEATURE_NUM;f2++){
            int missed = 0,missed_p=0,missed_n=0,g=0;//,correct=0;
            for(; g<sample_pCnt; g++){
                if(sample_f[f1][g]+sample_f[f2][g] <= 0){
                        missed_p++;
                }
            }
            for(; g<sample_allCnt;g++){
                if(sample_f[f1][g]+sample_f[f2][g] <= 0){
                        missed_n++;
                }
            }
            
            missed = missed_n*sample_w+missed_p*pos_weight;
            //if(WEIGHTED)
             //   missed = missed_n*(GENE_NUM-pos_num)/sample_allCnt+missed_p*pos_weight;
            //else
            //    missed = missed_n*GENE_NUM/sample_allCnt;
            if(missed+interval < best[0].miss){
                top_2d curFpair(f1,f2,missed+interval,missed_p,missed_n);
                best[0] = curFpair;
                sort(best,best+TOPK);
            }
            if(missed-interval < best[0].miss){
                top_2d curFpair(f1,f2,missed-interval,missed_p,missed_n);
                candidates.push_back(curFpair);
            }
            total_checked_g+=g;
            total_fp++;
        }
    }
    
    int per_checked_g = 2*total_checked_g/(2*FEATURE_NUM-FEATURE_CONSIDER-1)/(FEATURE_CONSIDER);
    cout << "total_checked_g="<<total_checked_g<<"; per_checked_g="<< per_checked_g << endl;
    
    fprintf(stderr,"------candidate size: %ld\n-------",candidates.size());
    fprintf(stderr,"------total feature pair considered: %d\n-------",total_fp);
    vector<top_2d > finalCandidates;
    for(int i=0; i<candidates.size();i++){
        if(best[0].miss > candidates[i].miss){
            finalCandidates.push_back(candidates[i]);
        } 
    }
    fprintf(stderr,"------finalCandidates size: %ld\n-------",finalCandidates.size());
    t = clock() - t;
    fprintf (stderr, "sampling_candidate_generation_horizontal: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    fprintf(ftime,"%ld\t%ld\t%d\t%f\t",finalCandidates.size(),total_checked_g,per_checked_g,((float)t)/CLOCKS_PER_SEC);

    return finalCandidates;
}

vector<top_2d > sampling_candidate_generation_vertical(FILE* ftime){
    clock_t t = clock();

    top_2d best[TOPK];
    vector<top_2d > candidates;
    for(int i=0;i<TOPK; i++){
        best[i].miss=2*GENE_NUM;
    }
    long total_checked_g=0;
    int sample_pCnt = (WEIGHTED==true)? pos_num:0; 
    int sample_nCnt = 1/(DELTA*DELTA);
    int sample_allCnt = sample_nCnt+sample_pCnt; 
    int interval = (WEIGHTED==true)?(GENE_NUM-pos_num)*DELTA:GENE_NUM*DELTA;  // confidence interval
    int sample_w = (WEIGHTED == true)? (GENE_NUM-pos_num)/sample_nCnt: GENE_NUM/sample_nCnt;

    //WEIGHTED: (GENE_NUM-pos_num)*DELTA; UNWEIGHTED:GENE_NUM*DELTA
    cout<<"sample_pCnt="<<sample_pCnt <<";sample_allCnt="<<sample_allCnt<<endl;
    cout << "confidence interval = " << interval << endl;
    int total_fp = 0;
    for(int f1=1;f1<FEATURE_CONSIDER;f1++){
        for(int f2=0; f2<f1;f2++){
            int missed = 0,missed_p=0,missed_n=0,g=0;//,correct=0;
            for(; g<sample_pCnt; g++){
                if(sample_f[f1][g]+sample_f[f2][g] <= 0){
                        missed_p++;
                }
            }
            for(; g<sample_allCnt;g++){
                if(sample_f[f1][g]+sample_f[f2][g] <= 0){
                        missed_n++;
                }
            }
            
            missed = missed_n*sample_w+missed_p*pos_weight;
            //if(WEIGHTED)
             //   missed = missed_n*(GENE_NUM-pos_num)/sample_allCnt+missed_p*pos_weight;
            //else
            //    missed = missed_n*GENE_NUM/sample_allCnt;
            if(missed+interval < best[0].miss){
                top_2d curFpair(f1,f2,missed+interval,missed_p,missed_n);
                best[0] = curFpair;
                sort(best,best+TOPK);
            }
            if(missed-interval < best[0].miss){
                top_2d curFpair(f1,f2,missed-interval,missed_p,missed_n);
                candidates.push_back(curFpair);
            }
            total_checked_g+=g;
            total_fp++;
        }
    }
    
    int per_checked_g = 2*total_checked_g/FEATURE_CONSIDER/(FEATURE_CONSIDER-1);
    cout << "total_checked_g="<<total_checked_g<<"; per_checked_g="<< per_checked_g << endl;
    
    fprintf(stderr,"------candidate size: %ld\n-------",candidates.size());
    fprintf(stderr,"------total feature pair considered: %d\n-------",total_fp);
    vector<top_2d > finalCandidates;
    for(int i=0; i<candidates.size();i++){
        if(best[0].miss > candidates[i].miss){
            finalCandidates.push_back(candidates[i]);
        } 
    }
    fprintf(stderr,"------finalCandidates size: %ld\n-------",finalCandidates.size());
    t = clock() - t;
    fprintf (stderr, "sampling_candidate_generation_vertical: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    fprintf(ftime,"%ld\t%ld\t%d\t%f\t",finalCandidates.size(),total_checked_g,per_checked_g,((float)t)/CLOCKS_PER_SEC);

    return finalCandidates;
}

top_2d calculate_real_error(int f1, int f2){
    int missed_p=0,missed_n=0,missed=0;
    for(int g=0; g<pos_num; g++){
        if(ff[f1].vals[g]+ff[f2].vals[g] <= 0)
            missed_p++;
    }
    for(int g=pos_num; g<GENE_NUM; g++){
        if(ff[f1].vals[g]+ff[f2].vals[g] <= 0)
            missed_n++;
    }
    missed = missed_p*pos_weight+missed_n;
    top_2d res(f1,f2,missed,missed_p,missed_n);
    return res;
    //fprintf(fresult,"\t%d\t%d\t%d\n",missed,missed_p,missed_n);
}

void two_phase_sampling(vector<top_2d > candidates, FILE* ftime){
    clock_t t = clock();
    top_2d real[TOPK];

    for(int i=0;i<TOPK; i++){
        real[i].miss=2*GENE_NUM;
    }

    for(int i=0; i<candidates.size();i++){
        top_2d cur = calculate_real_error(candidates[i].f1,candidates[i].f2);
        if(cur.miss < real[0].miss){
            real[0] = cur;
            sort(real,real+TOPK);
        }
    }    
    t = clock() - t;
    fprintf (stderr, "******best feature pair (%d,%d:%d) *****\n",ff[real[TOPK-1].f1].fid,ff[real[TOPK-1].f2].fid,real[TOPK-1].miss);
    fprintf (stderr, "two_phase_sampling find best feature pair: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
   
    FILE *real_file = fopen("FP_sampling.txt","a");
    for(int i=TOPK-1;i>=0;i--){
        int f1_tmp= ff[real[i].f1].fid, f2_tmp=ff[real[i].f2].fid;
        int f1 = min(f1_tmp,f2_tmp), f2 = max(f1_tmp,f2_tmp);
        fprintf(real_file,"%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", f1,f2,featureName[f1].c_str(),featureName[f2].c_str(),pos_num*pos_weight+(GENE_NUM-pos_num)-real[i].miss, real[i].miss,pos_num-real[i].miss_p, GENE_NUM-pos_num-real[i].miss_n, real[i].miss_p,real[i].miss_n);
    }
    fclose(real_file);
    fprintf(ftime,"%f\n",((float)t)/CLOCKS_PER_SEC);
}

void two_phase_sampling_opt(vector<top_2d > candidates, FILE* ftime){
    clock_t t = clock();
    top_2d real[TOPK];

    sort(candidates.begin(),candidates.end());
    cout << candidates[0].miss  <<" "<<candidates[100].miss <<" "<< candidates[candidates.size()-1].miss << endl;
    for(int i=0;i<TOPK; i++){
        real[i].miss=2*GENE_NUM;
    }

    int skip=0;
    for(int i=candidates.size()-1; i>=0;i--){
        if(real[0].miss < candidates[i].miss){
            skip++;
            continue;  // break;
        }
        top_2d cur = calculate_real_error(candidates[i].f1,candidates[i].f2);
        if(cur.miss < real[0].miss){
            real[0] = cur;
            sort(real,real+TOPK);
        }
    }    
    t = clock() - t;
    cout << "final candidates(opt)  = " << candidates.size()-skip << endl;
    fprintf (stderr, "******best feature pair (%d,%d:%d) *****\n",ff[real[TOPK-1].f1].fid,ff[real[TOPK-1].f2].fid,real[TOPK-1].miss);
    fprintf (stderr, "two_phase_sampling_opt find best feature pair: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
   
    FILE *real_file = fopen("FP_sampling_opt.txt","a");
    for(int i=TOPK-1;i>=0;i--){
        int f1_tmp= ff[real[i].f1].fid, f2_tmp=ff[real[i].f2].fid;
        int f1 = min(f1_tmp,f2_tmp), f2 = max(f1_tmp,f2_tmp);
        fprintf(real_file,"%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", f1,f2,featureName[f1].c_str(),featureName[f2].c_str(),pos_num*pos_weight+(GENE_NUM-pos_num)-real[i].miss, real[i].miss,pos_num-real[i].miss_p, GENE_NUM-pos_num-real[i].miss_n, real[i].miss_p,real[i].miss_n);        
    }
    fclose(real_file);
    fprintf(ftime,"%ld\t%f\n",candidates.size()-skip,((float)t)/CLOCKS_PER_SEC);
}


void print_fp(int f1, int f2, int cnt){
    clock_t t = clock();
   
    /*int missed = 0,missed_p=0,missed_n=0;//,correct=0;
    int g=0;
    for(; g<pos_num; g++){
        if(ff[f1].vals[g]+ff[f2].vals[g] <= 0){
            missed_p++;
        }
    }
    for(; g<GENE_NUM; g++){
        if(ff[f1].vals[g]+ff[f2].vals[g] <= 0){    
            missed_n++;
        }
    }
    missed = missed_p*pos_weight+missed_n;
    t = clock() - t;
    
    fprintf (stderr, "******%d,%d(%d,%d:%d)*****\n",f1,f2,missed_p,missed_n,missed);
    FILE *analysis_file = fopen("FPanalysis.txt","w");
    
    for(int g=0; g<MATRIX_GENE_NUM;g++){
        if(genes[g]==1)
            fprintf(analysis_file,"%f\t%f\n",features[f1].vals[g],features[f2].vals[g]);
    }
    fprintf(analysis_file,"\n\n");
    for(int g=0; g<MATRIX_GENE_NUM;g++){
        if(genes[g]==-1)
            fprintf(analysis_file,"%f\t%f\n",features[f1].vals[g],features[f2].vals[g]);
    }
    fclose(analysis_file);*/
    
    val_type xmin, xmax, ymin, ymax;
    vector<point > pos;
    vector<point > neg;
    vector<val_type > xx = features[f1].vals;
    vector<val_type > yy = features[f2].vals;
    for(int i=0; i<pos_num; i++){
        point cur(features[f1].vals[i],features[f2].vals[i]);
        pos.push_back(cur);
    }
    for(int i=pos_num; i<GENE_NUM; i++){
        point cur(features[f1].vals[i],features[f2].vals[i]);
        neg.push_back(cur);
    }
    
    sort(xx.begin(),xx.end());
    sort(yy.begin(),yy.end());
    xmax = xx[0.95*GENE_NUM]; //GENE_NUM-1];  // 0.95*GENE_NUM
    xmin = xx[0.05*GENE_NUM];
    ymin = yy[0.05*GENE_NUM];//0];
   	ymax = yy[0.95*GENE_NUM];
   	
   	double xrange = xmax- xmin;
   	double yrange = ymax- ymin;
    
    cout << xmin << " " << xmax << " " << ymin <<" " << ymax << endl;
    
    vector<vector<int> > hist_pos(cnt,vector<int>(cnt,0));
    for(int i=0; i<pos.size();i++){
    	int ptrx = (pos[i].x-xmin)*cnt/xrange;
    	int ptry = (pos[i].y-ymin)*cnt/yrange;
    	if(ptrx>=cnt) ptrx = cnt-1;
    	if(ptry>=cnt) ptry = cnt-1;
    	if(ptrx<0) ptrx = 0;
    	if(ptry<0) ptry = 0;
    	hist_pos[ptrx][ptry]+=pos_weight;
    }
    FILE *fpos = fopen("posHist.txt","w");
    for(int i=0; i<cnt; i++){
	    fprintf(fpos, "\t%lf", xmin+xrange*i/cnt);    	
    }
    fprintf(fpos, "\n");
    int test_pos=0;
    for(int i=0; i<cnt;i++){
    	fprintf(fpos, "%lf", ymin+yrange*i/cnt);
    	for(int j=0; j<cnt; j++){
    		fprintf(fpos, "\t%d", hist_pos[i][j]);
    		test_pos+=hist_pos[i][j];
    	}
    	fprintf(fpos, "\n");
    }
    fclose(fpos);
    cout << "test_pos" << test_pos << endl;
     
    vector<vector<int> > hist_neg(cnt,vector<int>(cnt,0));
    for(int i=0; i<neg.size();i++){
    	int ptrx = (neg[i].x-xmin)*cnt/xrange;
    	int ptry = (neg[i].y-ymin)*cnt/yrange;
    	if(ptrx>=cnt) ptrx = cnt-1;
    	if(ptry>=cnt) ptry = cnt-1;
    	if(ptrx<0) ptrx = 0;
    	if(ptry<0) ptry = 0;
    	if(ptrx >cnt || ptry >cnt || ptry <0 || ptrx <0)
    		cout << "error======" << ptrx<< " " << ptry <<" " << neg[i].x <<" " <<neg[i].y <<endl;
    	hist_neg[ptrx][ptry]++;
    }
    
    FILE *fneg = fopen("negHist.txt","w");
	for(int i=0; i<cnt; i++){
	    fprintf(fneg, "\t%lf", xmin+xrange*i/cnt);    	
    }
    fprintf(fneg, "\n");

    int test_neg=0;
    for(int i=0; i<cnt;i++){
    	fprintf(fneg, "%lf", ymin+yrange*i/cnt);
    	for(int j=0; j<cnt; j++){
    		fprintf(fneg, "\t%d", hist_neg[i][j]);
    		test_neg+=hist_neg[i][j];
    	}
    	fprintf(fneg, "\n");
    }
    fclose(fneg);
    cout << "test_neg" << test_neg << endl;

    FILE *fratio = fopen("ratioHist.txt","w");
    for(int i=0; i<cnt; i++){
        fprintf(fratio, "\t%lf", xmin+xrange*i/cnt);      
    }
    fprintf(fratio, "\n");
    for(int i=0; i<cnt;i++){
        fprintf(fratio, "%lf", ymin+yrange*i/cnt);
        for(int j=0; j<cnt; j++){
            //double ratio = (double)(hist_pos[i][j]-hist_neg[i][j])/(hist_pos[i][j]+hist_neg[i][j]);
            //fprintf(fratio, "\t%f", ratio);
            int diff = hist_pos[i][j]-hist_neg[i][j];
            int diffsqr = (diff==0?1:sqrt(abs(diff)));
            fprintf(fratio, "\t%d", diff/diffsqr);
        }
        fprintf(fratio, "\n");
    }
    fclose(fratio);
    t = clock() - t;
    fprintf (stderr, "It took me %ld clicks (%f seconds) for the histgram program.\n",t,((float)t)/CLOCKS_PER_SEC);
}

int main(int argc, char* argv[]){
    FILE *ftime=fopen("Time.txt","a");
    initlization();
    parse_args(argc, argv);
    cout << "SORT_F=" << SORT_F << endl;
    
    if(HIST){
        if(DELIMITER.compare("TRANSPOSE") == 0 ){
            histogram_sample_feature_matrix_space();
        }else{
            if(DELIMITER.compare("COMMA") == 0)
                histogram_feature_sample_matrix_comma();
            else
                histogram_feature_sample_matrix_space();
        }
        return 0;
    }
    
    if(DELIMITER.compare("COMMA") == 0 )
        load_feature_gene_matrix_comma();
    else{
        if(DELIMITER.compare("TRANSPOSE") == 0 )
            load_sample_feature_matrix_space();
        else
            load_feature_gene_matrix_space();
    }
        
    if(!REST)  // P v.s. N or p+n v.s. rest
        load_exp_pos_neg();  // needs change!!!
    else
        load_exp_rest();
    pos_weight=1;
    if(WEIGHTED == true)
        pos_weight = (GENE_NUM-pos_num)/pos_num;
    cout << "pos_weight=" << pos_weight << endl;
    if(FEATURE_CONSIDER_TOTAL == 0)
        FEATURE_CONSIDER = FEATURE_NUM;
    else{
        if(HORIZONTAL)
            FEATURE_CONSIDER = ((2*FEATURE_NUM-1)-sqrt(((long)(2*FEATURE_NUM-1))*(2*FEATURE_NUM-1)-8*FEATURE_CONSIDER_TOTAL))/2;
        else
            FEATURE_CONSIDER = (1+sqrt(1+8*FEATURE_CONSIDER_TOTAL))/2;
    }
    cout << "FEATURE_NUM=" << FEATURE_NUM << "FEATURE_CONSIDER=" <<FEATURE_CONSIDER << "; GENE_NUM =" << GENE_NUM << endl;
    
    if(!TRASNFORM){
        cout << "no transformation" << endl;
        no_transformation_baseline(ftime);
    }else{
        if(PRINT){
            cout << "printF1=" <<printF1 <<";printF2=" << printF2 << endl;
            print_fp(printF1,printF2,buckets);
        }else{
            transformation();
            //fprintf(stderr,"-----done with transformation------\n");
            cout << "FEATURE_CONSIDER=" << FEATURE_CONSIDER <<endl;
            if(SORT_G==true){
                if(WEIGHTED==true)
                    sort_genes_weighted();
                else
                    sort_genes();
            }
            if(SORT_F==true){
                sort(ff.begin(),ff.end());
                cout <<"***" << ff[0].correct <<" " << ff[1].correct << " " << ff[ff.size()-1].correct << "***" << endl;
            }
            if(EARLY_TERMINATE && HORIZONTAL)
                early_term_horizontal(ftime);
            if(EARLY_TERMINATE && !HORIZONTAL)
                early_term_vertical(ftime);
            
            if(SAMPLING || SAMPLING_OPT){
                //fprintf(stderr,"-----here------\n");
                prepare_samples();
                vector<top_2d > candidates;
                if(HORIZONTAL)
                    candidates = sampling_candidate_generation_horizontal(ftime);
                else
                    candidates = sampling_candidate_generation_vertical(ftime);
                if(SAMPLING)
                    two_phase_sampling(candidates,ftime);   
                else
                    two_phase_sampling_opt(candidates,ftime);
            }else
            if(!EARLY_TERMINATE){  // no sortG no sortF
                if(HORIZONTAL)
                    baseline(ftime);
                else
                    baseline_vertical(ftime);
            }
        }
    }
}
