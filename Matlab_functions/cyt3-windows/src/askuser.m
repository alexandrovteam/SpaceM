function y=askuser(s)
    choice = questdlg(s, 'Question', 'Yes', 'No', 'Cancel', 'Cancel');
    y = strcmpi(choice, 'Yes');
end