function[table2Save] = saveCoefficientsAndErrors(coeffss, coeffNamess, errs)



%%% Put the coefficient names into a table

valuesAndErrors = zeros(size(coeffss,2),2);

for i =1:size(coeffss,2)
  

   valuesAndErrors(i,1) = coeffss(i);
   valuesAndErrors(i,2) = errs(i); 
end
table2Save = array2table(valuesAndErrors,'VariableNames',{'Values','Errors'},...
    'RowNames',coeffNamess);

%writetable(table2Save,fileNameSave,'WriteRowNames',true);


end

