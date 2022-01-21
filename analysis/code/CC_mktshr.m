function [prod_shr, firm_shr] = CC_mktshr(Z)
    % Take Z computed from CC_WIODaug.m and generate firm level market
    % share data across HF, PR, and chip.
    
    % Market shares in each country.
    Z_HF    = Z(:,:,end,end-2)./sum(Z(:,:,end,end-2));
    Z_PR    = Z(:,:,end-1,end-2)./sum(Z(:,:,end-1,end-2));
    Z_chip  = Z(:,:,end-2,end-3)./sum(Z(:,:,end-2,end-3));
    
    % First, assume HF markets are divided into five producers
    n =5;
    
    prod_shr.HFshr = kron(Z_HF,ones(n,1)/5); 
    firm_shr.HFshr = ones(size(Z_HF,1)*n,size(Z_HF,2))/n;
    % Next, take the actual market share data of PR.
    % 2014 KrF PR global market share
    PR_shr = [28.2 20.0 17.3 14.5 8.2 6.4 5.4]'/100;
    prod_shr.PRshr = kron(Z_PR,PR_shr);
    firm_shr.PRshr = repmat(PR_shr,size(Z_PR,1),size(Z_PR,2));
    firm_shr.PRshr(prod_shr.PRshr==0)=0;
    % Lastly, compute the chip market share and compute S(j)^(k).
    Dram_shr = [40.4 27.4 24.6 0 0 0]'; % global D-ram market share (IHS)
    % Samsung, SK hynix, and Micron.
    % Dram market size 2014: 46.9 billion (IDC)
    % NAND flash market size 2014: 28.9 billion
    NAND_shr = [30.7 11.8 14.6 19.3 15.6 6.8]'; %source: (IHS)
    % Samsung, SK hynix, Micron, Western digital, Kioxia, Intel 
    
    Dram_sales = (Dram_shr/100)*46.9;
    NAND_sales = (NAND_shr/100)*28.9;
    chip_sales = Dram_sales+NAND_sales;
    chip_shr = chip_sales/sum(chip_sales);
    prod_shr.chipshr = kron(Z_chip,chip_shr);
    firm_shr.chipshr = repmat(chip_shr,size(Z_chip,1),size(Z_chip,2));
end