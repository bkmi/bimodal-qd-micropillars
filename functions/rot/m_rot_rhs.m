function [ y ] = m_rot_rhs(xx,p,A,expA,user_rhs,user_tau,isvec)
%My own version of the right hand side eqn.
%   Input:
%       xx, ...
%       p, ...
%       A, ...
%       expA, ...
%       user_rhs, ...
%       user_tau, ...
%       isvec
%   Output:
%       y

%% right-hand side in rotating coordinates
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%

if isa(A,'double')
    % USED TO WORK THIS WAY, NEW VERSION BELOW IS FOR MULTIVARIABLE EXPA
    omega=p(end);
    userpar=p(1:end-1);
    xxrot=xx;
    tau_ind=user_tau();
    if ~isvec
        for i=2:size(xx,2)
            xxrot(:,i)=expA(-omega*p(tau_ind(i-1)))*xxrot(:,i);
        end
    else
        dim=size(xxrot,1); % Number of rows
        nvec=size(xxrot,3); % Number of 'sheets' -- vectorized
        for i=2:size(xx,2) % 2:(number of columns)
    %         xxrot(:,i,:) % all rows, second column, all pages 
    %         expA(-omega*p(tau_ind(i-1))) %expA @ tau_fb1
    %         reshape(xxrot(:,i,:),dim,nvec) % Turn multiple sheets into a 2d matrix
            xxrot(:,i,:)=expA(-omega*p(tau_ind(i-1)))*reshape(xxrot(:,i,:),dim,nvec);
            % multiply the matricies and move the resulting columns into the
            % relevant page.
            % 
            % The whole process extracted the second column from each page, did
            % the multiplication, then put each multiplied second column back
            % where it came from.
        end
    end
    y0=user_rhs(xxrot,userpar); % calculate: given new, rotated 2nd column
    % y0 == f(y(t),exp(-A*omega*tau1)*y(t-tau1),...)


    % A*omega*reshape(xxrot(:,1,:),[dim,nvec])
    % This thing is the -A*w*y(t) part

    % reshape(A*omega*reshape(xxrot(:,1,:),[dim,nvec]),[dim,1,nvec])
    % This thing puts it back into pages

    y=y0-reshape(A*omega*reshape(xxrot(:,1,:),[dim,nvec]),[dim,1,nvec]);
elseif isa(A,'cell')
    % Arbitrary rotations and expA(_,_,...,_) s.t. expA has numel(A)
    % inputs. Aka. expA takes in a number of parameters equal to the number
    % of rotation matricies you use.
    
    numRot = numel(A); % Number of rot params/matricies.
    omega = cell(numRot,1);
    for i = 1:numRot
        omega{i} = p(end-numRot+i);
    end
    userpar = p(1:end-numRot);
    xxrot = xx;
    tau_ind = user_tau();

    if ~isvec
        for i=2:size(xx,2)
            % Updates the ith column of xxrot evaluated at the relevant
            % delay. i starts at 2 (the first delay).
            
            % Now I populate a new cell which contains each omega
            % multiplied by the relevant delay.
            omegaMulDelay = num2cell( ... % Convert to cell for later.
                cellfun(@(x)-x*p(tau_ind(i-1)),omega) ... % Each Omega * tau(i-1)
                );
            
            xxrot(:,i) = ...
                expA( omegaMulDelay{:} ) ...
                * xxrot(:,i);
        end
    else
        dim=size(xxrot,1); % Number of rows
        nvec=size(xxrot,3); % Number of 'sheets' -- vectorized
        for i=2:size(xx,2) % 2:(number of columns)
            % Updates the ith column of xxrot evaluated at the relevant
            % delay. i starts at 2 (the first delay).
            
            % Now I populate a new cell which contains each omega
            % multiplied by the relevant delay.
            omegaMulDelay = num2cell( ... % Convert to cell for later.
                cellfun(@(x)-x*p(tau_ind(i-1)),omega) ... % Each Omega * tau(i-1)
                );
            xxrot(:,i,:) = ...
                expA(omegaMulDelay{:}) ... 
                *reshape(xxrot(:,i,:),dim,nvec);
            % multiply the matricies and move the resulting columns into
            % the relevant page.
            % 
            % The whole process extracted the second column from each page,
            % did the multiplication, then put each multiplied second
            % column back where it came from... It's just organizational.
        end
    end
    y0=user_rhs(xxrot,userpar); % calculate: given new, rotated 2nd column
    % y0 == f(y(t),exp(-A*omega*tau1)*y(t-tau1),...) <-- the math
    
    % A*omega*reshape(xxrot(:,1,:),[dim,nvec])
    % This thing is the -A*w*y(t) part

    % reshape(A*omega*reshape(xxrot(:,1,:),[dim,nvec]),[dim,1,nvec])
    % This thing puts it back into pages
    
    % We need to multiple the ith omega with the ith rotation. I will do
    % this first in a loop.
    omegaMulRotationSum = 0;
    for i = 1:numRot
        omegaMulRotationSum = omegaMulRotationSum + A{i}*omega{i};
    end
    
    y = y0 ...
        - reshape( ...
        ( omegaMulRotationSum ) ... % Sum (Rotation matrix * relevant Omega)
        * reshape(xxrot(:,1,:),[dim,nvec]), ... % Second reshape
        [dim,1,nvec] ...
        ); % Finish first reshape
    
%%    GIVEN TWO ROTATIONS AND EXPA WITH TWO VAR
%     
%     % userpar, omega1, omega2, omega3, ...
%     omega1=p(end-1);
%     omega2=p(end);
%     userpar=p(1:end-2);
%     xxrot=xx;
%     tau_ind=user_tau();
%     
%     if ~isvec
%         for i=2:size(xx,2)
%             xxrot(:,i)=expA(-omega1*p(tau_ind(i-1)),-omega2*p(tau_ind(i-1)))*xxrot(:,i);
%         end
%     else
%         dim=size(xxrot,1); % Number of rows
%         nvec=size(xxrot,3); % Number of 'sheets' -- vectorized
%         for i=2:size(xx,2) % 2:(number of columns)
% %             xxrot(:,i,:) % all rows, ith column, all pages
% %             expA(-omega*p(tau_ind(i-1))) %expA @ tau_fb1 (One omega)
% %             reshape(xxrot(:,i,:),dim,nvec) % Turn multiple sheets into a
% %             2d matrix
%             xxrot(:,i,:)=expA(-omega1*p(tau_ind(i-1)),-omega2*p(tau_ind(i-1)))*reshape(xxrot(:,i,:),dim,nvec);
%             % multiply the matricies and move the resulting columns into
%             % the relevant page.
%             % 
%             % The whole process extracted the second column from each page,
%             % did the multiplication, then put each multiplied second
%             % column back where it came from.
%         end
%     end
%     y0=user_rhs(xxrot,userpar); % calculate: given new, rotated 2nd column
%     % y0 == f(y(t),exp(-A*omega*tau1)*y(t-tau1),...)
% 
% 
%     % A*omega*reshape(xxrot(:,1,:),[dim,nvec])
%     % This thing is the -A*w*y(t) part
% 
%     % reshape(A*omega*reshape(xxrot(:,1,:),[dim,nvec]),[dim,1,nvec])
%     % This thing puts it back into pages
% 
%     y=y0-reshape((A{1}*omega1+A{2}*omega2)*reshape(xxrot(:,1,:),[dim,nvec]),[dim,1,nvec]);
end




end

