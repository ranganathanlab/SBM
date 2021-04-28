function [h_out,s_out,y_out,ys_out,hess_diag_out,ind_out,skipping] = UpdateHessian(g,y,s,s_out,y_out,ys_out,hess_diag,options,ind,skipping)

        ys = y'*s;  
        if ys > 1e-10
            y_out(:,(ind==max(ind)))=y;

            s_out(:,(ind==max(ind)))=s;
            ys_out((ind==max(ind)))=ys;

            hess_diag_out = ys/(y'*y);
            
            if ind(options.m)==0
                ind_out=ind;
                ind_out(find(ind==max(ind))+1)=max(ind)+1;
            else 
                ind_out=ind([end 1:end-1]);
            end
                
            
        else

            hess_diag_out=hess_diag;
            ind_out=ind;
            skipping=skipping+1;
        end
        
        [~,order]=sort(ind(ind>0));
        
        % calculate the -gradient*inv_hessian
        h_out = -g;
        for i=(order(end:-1:1))
            alpha(i) = (s_out(:,i)'*h_out)/ys_out(i);
            h_out = h_out-alpha(i)*y_out(:,i);
        end
        h_out = hess_diag_out*h_out;
        for i=order
            beta(i) = (y_out(:,i)'*h_out)/ys_out(i);
            h_out = h_out + s_out(:,i)*(alpha(i)-beta(i));
        end
end

