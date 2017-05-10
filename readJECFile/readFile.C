#include <string>
#include <fstream>
#include <iostream>

std::string readFile(const char *filename = "80X_mcRun2_asymptotic_2016_miniAODv2_v1_L2Relative_AK4Calo.txt", int linenumber = 42) {

  std::ifstream infile(filename) ;

  int lines_read = 0 ;
  
  std::string line ;
   if ( infile.is_open( ) ) {
      while ( infile ) {
	 getline( infile , line ) ;
	 lines_read++ ;
	 if ( lines_read == linenumber ) {
	    std::cout << line << std::endl ;
            // float par[12];
            // sscanf(line.c_str(),"%f %f %f %f %f %f %f %f %f %f %f %f",&par[0],&par[1],&par[2],&par[3],&par[4],&par[5],&par[6],&par[7],&par[8],&par[9],&par[10],&par[11]);
            // for(int i = 0; i<12; ++i)
            //   Printf("par[%d] = %f",i,par[i]);
	    break ; 
	 }
      }
      infile.close( ) ;
      if ( lines_read < linenumber ) 
	 std::cout << "No " << linenumber << " lines in " << filename << " !\n" ;
      return line ;
   }
   else {
      std::cerr << "Could not find file " << filename << " !\n" ;
      return "" ;
   }

}
