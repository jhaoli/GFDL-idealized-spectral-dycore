  �$  V   k820309    �          17.0        a _                                                                                                           
       /public/home/wab/lijh/GFDL/full_dynamics/../src/atmos_spectral/model/matrix_invert.F90 MATRIX_INVERT_MOD              MAXMAG                      @                              
       MPP_PE MPP_ROOT_PE ERROR_MESG FATAL WRITE_VERSION_NUMBER                �  @                              '                    #BEGIN    #END    #SIZE    #MAX_SIZE    #IS_GLOBAL                 � D                                                             � D                                                            � D                                                            � D                                                            � D                                                              �  @                              '�                   #NAME 	   #UNITS 
   #LONGNAME    #STANDARD_NAME    #MIN    #MAX    #MISSING    #FILL    #SCALE    #ADD    #PACK    #AXES    #SIZE 6   #TIME_AXIS_INDEX 7   #ID 8   #TYPE 9   #NATT :   #NDIM ;   #ATT <   #POSITION =                � D                            	     �                                   � D                            
     �       �                           � D                                                                   � D                                 �                                  � D                                  �         
                � D                                  �         
                � D                                  �         
                � D                                  �         
                � D                                  �      	   
                � D                                  �      
   
                � D                                  �                      �D                                         �      �           #AXISTYPE              &                                                                   P              y#AXISTYPE                                                                  �  @                              '�                   #NAME    #UNITS    #LONGNAME    #CARTESIAN    #CALENDAR    #SENSE    #LEN    #DOMAIN    #DATA )   #ID *   #DID +   #TYPE ,   #NATT -   #SHIFT .   #ATT /                � D                                 �                                   � D                                 �       �                           � D                                                                   � D                                                                   � D                                                                  � D                                                            � D                                  $                         � D                                  �       (             #DOMAIN1D                   �  @                              '�              
      #COMPUTE    #DATA     #GLOBAL !   #MEMORY "   #CYCLIC #   #LIST $   #PE %   #POS &   #GOFFSET '   #LOFFSET (                � D                                                        #DOMAIN_AXIS_SPEC                 � D                                                        #DOMAIN_AXIS_SPEC                 � D                             !            (              #DOMAIN_AXIS_SPEC                 � D                             "            <              #DOMAIN_AXIS_SPEC                 � D                             #     P                       �D                             $            X       �            #DOMAIN1D              &                                                                   P              y#DOMAIN1D                                                                � D                             %     �                          � D                             &     �                          � D                             '     �       	                   � D                             (     �       
                �D                             )            �             	  
            &                                                                   P              y
                                                            � D                             *            
                   � D                             +     $                         � D                             ,     (                         � D                             -     ,                         � D                             .     0                      �D                             /            8      �           #ATTTYPE 0             &                                                                   P              y#ATTTYPE 0                                                                 �  @                         0     '�                   #TYPE 1   #LEN 2   #NAME 3   #CATT 4   #FATT 5                � D                             1                                � D                             2                               � D                            3     �                                  � D                            4            �                        �D                             5            �               
            &                                                                   P              y
                                                          �D                             6                                        &                                                                   P              y                                                            � D                             7     H                         � D                             8     L                         � D                             9     P                         � D                             :     T                         � D                             ;     X                      �D                             <            `      �           #ATTTYPE 0             &                                                                   P              y#ATTTYPE 0                                                               � D                             =     �            %         @                                 >                            %         @                                 ?                            #         @                                  @                    #ROUTINE A   #MESSAGE B   #LEVEL C             
                                 A                    1           
                                 B                    1           
                                  C                                                        D                                                      2#         @                                  E                    #VERSION F   #TAG G   #UNIT H             
                                 F                    1           
                                G                    1           
                                 H                                                   I                   
                &                                                                                   J                   
                &                                           #         @                                   K                    #MATRIX L   #DET M          0  
D@                               L                   
               &                   &                                                     D @                               M     
                     @                              N                                                        O     �                            �                       C$Id matrix_invert.f90 $                                                                                                                                                                     P     �                            �                       C$Name: siena_201207 $                                                                                                                                     @                                 Q            %         @                                R                           #H S   #M T            
                                  S                    
    p          5 � p        r T       5 � p        r T                               
                                  T              �   q      fn#fn '        b   uapp(MATRIX_INVERT_MOD    (  y   J  FMS_MOD 1   �  �       DOMAIN_AXIS_SPEC+MPP_DOMAINS_MOD =   ,  H   %   DOMAIN_AXIS_SPEC%BEGIN+MPP_DOMAINS_MOD=BEGIN 9   t  H   %   DOMAIN_AXIS_SPEC%END+MPP_DOMAINS_MOD=END ;   �  H   %   DOMAIN_AXIS_SPEC%SIZE+MPP_DOMAINS_MOD=SIZE C     H   %   DOMAIN_AXIS_SPEC%MAX_SIZE+MPP_DOMAINS_MOD=MAX_SIZE E   L  H   %   DOMAIN_AXIS_SPEC%IS_GLOBAL+MPP_DOMAINS_MOD=IS_GLOBAL %   �  3      FIELDTYPE+MPP_IO_MOD /   �  P   %   FIELDTYPE%NAME+MPP_IO_MOD=NAME 1     P   %   FIELDTYPE%UNITS+MPP_IO_MOD=UNITS 7   g  P   %   FIELDTYPE%LONGNAME+MPP_IO_MOD=LONGNAME A   �  P   %   FIELDTYPE%STANDARD_NAME+MPP_IO_MOD=STANDARD_NAME -     H   %   FIELDTYPE%MIN+MPP_IO_MOD=MIN -   O  H   %   FIELDTYPE%MAX+MPP_IO_MOD=MAX 5   �  H   %   FIELDTYPE%MISSING+MPP_IO_MOD=MISSING /   �  H   %   FIELDTYPE%FILL+MPP_IO_MOD=FILL 1   '  H   %   FIELDTYPE%SCALE+MPP_IO_MOD=SCALE -   o  H   %   FIELDTYPE%ADD+MPP_IO_MOD=ADD /   �  H   %   FIELDTYPE%PACK+MPP_IO_MOD=PACK /   �    %   FIELDTYPE%AXES+MPP_IO_MOD=AXES $   	  �       AXISTYPE+MPP_IO_MOD .   
  P   %   AXISTYPE%NAME+MPP_IO_MOD=NAME 0   R
  P   %   AXISTYPE%UNITS+MPP_IO_MOD=UNITS 6   �
  P   %   AXISTYPE%LONGNAME+MPP_IO_MOD=LONGNAME 8   �
  P   %   AXISTYPE%CARTESIAN+MPP_IO_MOD=CARTESIAN 6   B  P   %   AXISTYPE%CALENDAR+MPP_IO_MOD=CALENDAR 0   �  H   %   AXISTYPE%SENSE+MPP_IO_MOD=SENSE ,   �  H   %   AXISTYPE%LEN+MPP_IO_MOD=LEN 2   "  ^   %   AXISTYPE%DOMAIN+MPP_IO_MOD=DOMAIN )   �  �       DOMAIN1D+MPP_DOMAINS_MOD 9   @  f   %   DOMAIN1D%COMPUTE+MPP_DOMAINS_MOD=COMPUTE 3   �  f   %   DOMAIN1D%DATA+MPP_DOMAINS_MOD=DATA 7     f   %   DOMAIN1D%GLOBAL+MPP_DOMAINS_MOD=GLOBAL 7   r  f   %   DOMAIN1D%MEMORY+MPP_DOMAINS_MOD=MEMORY 7   �  H   %   DOMAIN1D%CYCLIC+MPP_DOMAINS_MOD=CYCLIC 3        %   DOMAIN1D%LIST+MPP_DOMAINS_MOD=LIST /   0  H   %   DOMAIN1D%PE+MPP_DOMAINS_MOD=PE 1   x  H   %   DOMAIN1D%POS+MPP_DOMAINS_MOD=POS 9   �  H   %   DOMAIN1D%GOFFSET+MPP_DOMAINS_MOD=GOFFSET 9     H   %   DOMAIN1D%LOFFSET+MPP_DOMAINS_MOD=LOFFSET .   P  �   %   AXISTYPE%DATA+MPP_IO_MOD=DATA *   D  H   %   AXISTYPE%ID+MPP_IO_MOD=ID ,   �  H   %   AXISTYPE%DID+MPP_IO_MOD=DID .   �  H   %   AXISTYPE%TYPE+MPP_IO_MOD=TYPE .     H   %   AXISTYPE%NATT+MPP_IO_MOD=NATT 0   d  H   %   AXISTYPE%SHIFT+MPP_IO_MOD=SHIFT ,   �    %   AXISTYPE%ATT+MPP_IO_MOD=ATT #   �  �       ATTTYPE+MPP_IO_MOD -   ;  H   %   ATTTYPE%TYPE+MPP_IO_MOD=TYPE +   �  H   %   ATTTYPE%LEN+MPP_IO_MOD=LEN -   �  P   %   ATTTYPE%NAME+MPP_IO_MOD=NAME -     P   %   ATTTYPE%CATT+MPP_IO_MOD=CATT -   k  �   %   ATTTYPE%FATT+MPP_IO_MOD=FATT /   _  �   %   FIELDTYPE%SIZE+MPP_IO_MOD=SIZE E   S  H   %   FIELDTYPE%TIME_AXIS_INDEX+MPP_IO_MOD=TIME_AXIS_INDEX +   �  H   %   FIELDTYPE%ID+MPP_IO_MOD=ID /   �  H   %   FIELDTYPE%TYPE+MPP_IO_MOD=TYPE /   +  H   %   FIELDTYPE%NATT+MPP_IO_MOD=NATT /   s  H   %   FIELDTYPE%NDIM+MPP_IO_MOD=NDIM -   �    %   FIELDTYPE%ATT+MPP_IO_MOD=ATT 7   �  H   %   FIELDTYPE%POSITION+MPP_IO_MOD=POSITION      P       MPP_PE+MPP_MOD $   a  P       MPP_ROOT_PE+MPP_MOD #   �  m       ERROR_MESG+FMS_MOD +     L   a   ERROR_MESG%ROUTINE+FMS_MOD +   j  L   a   ERROR_MESG%MESSAGE+FMS_MOD )   �  @   a   ERROR_MESG%LEVEL+FMS_MOD (   �  q       FATAL+MPP_PARAMETER_MOD -   g  h       WRITE_VERSION_NUMBER+FMS_MOD 5   �  L   a   WRITE_VERSION_NUMBER%VERSION+FMS_MOD 1     L   a   WRITE_VERSION_NUMBER%TAG+FMS_MOD 2   g  @   a   WRITE_VERSION_NUMBER%UNIT+FMS_MOD /   �  �       MPP_DOMAINS_STACK+MPP_DATA_MOD 8   3  �       MPP_DOMAINS_STACK_NONBLOCK+MPP_DATA_MOD    �  ]       INVERT       �   a   INVERT%MATRIX    �   @   a   INVERT%DET     !  @       MAXMAG    @!        VERSION    A"        TAGNAME &   B#  @       ENTRY_TO_LOGFILE_DONE    �#  ^       MAX_MAG    �#  �   a   MAX_MAG%H    �$  @   a   MAX_MAG%M 