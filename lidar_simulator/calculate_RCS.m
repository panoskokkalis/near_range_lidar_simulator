function [RCS]=calculate_RCS(alt, elastic_data)

RCS = ones(size(elastic_data));

for i=1:1:size(elastic_data,2)
    RCS(:,i) = elastic_data(:,i).*(alt(:,1).^2);
end
