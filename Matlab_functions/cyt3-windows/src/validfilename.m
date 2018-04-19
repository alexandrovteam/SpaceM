function filename=validfilename(s)
  legalchars = 'a-zA-Z0-9\-\ \_\.' ;
  filename = regexprep(s,['[^' legalchars ']'],'_');
end