function dummy  = get_Gss()
DOCID = '1QFpSVeAzOmB7qCJ7PvhwqzQD9IQARc2d5su7NrHodu4';
ss = GetGoogleSpreadsheet(DOCID);
C = ss(2:end,1:5);
T = cell2table(C);
T.Properties.VariableNames = ss(1,1:5);
writetable(T, 'C:\Users\Luca\Documents\python_codebase\sf2an\Gss.csv');

dummy=0;
end
%      Given s:
%
%          s.Alpha = { 'First', 'Second';
%                      'Third', 'Fourth'};
%
%          s.Beta  = [[      1,       2;
%                            3,       4]];
%          
%          s.Gamma = {       1,       2;
%                            3,       4};
%
%          s.Epsln = [     abc;
%                          def;
%                          ghi];
% 
%      STRUCT2CSV(s,'any.csv') will produce a file 'any.csv' containing:
%
%         "Alpha",        , "Beta",   ,"Gamma",   , "Epsln",
%         "First","Second",      1,  2,      1,  2,   "abc",
%         "Third","Fourth",      3,  4,      3,  4,   "def",
%                ,        ,       ,   ,       ,   ,   "ghi",