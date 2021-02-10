function varargout = subsref(this,s)
switch s(1).type
    case '.'
        if length(s) == 1
            % Implement this.PropertyName
            varargout{1} = this.(s(1).subs);
            ...
        elseif length(s) == 2 && strcmp(s(2).type,'()')
            varargout{1} = this.(s(1).subs)(s(2).subs{1});
        else
        [varargout{1:nargout}] = builtin('subsref',this,s);
        end
    case '()'
        if length(s) == 1
            mc = metaclass(this) ; 
            indConstant = [mc.PropertyList.Constant] ; 
            indDependent = [mc.PropertyList.Dependent] ; 
            all_fields = {mc.PropertyList.Name} ; 
            fields = all_fields(~indConstant & ~indDependent) ; 
            for n=1:length(fields)
                if strcmp(fields{n}, 'gamma') 
                    tmpStruct.gamma = this.(fields{n}) ; 
                elseif strcmp(fields{n}, 'Nacq')
                    tmpStruct.Nacq = sum(s.subs{:}) ; 
                elseif isempty(this.(fields{n}))
                    tmpStruct.(fields{n}) = [] ; 
                else
                    tmpStruct.(fields{n}) = this.(fields{n})(s.subs{:},:) ; 
                end
            end
            % ----------------------- output a new object by calling a constractor ---------------------------
            % Note: for convenience, Please add all new subclass constractors here
            if isa(this,'mig.DiffusionPulseSequence')
                varargout{1} = mig.DiffusionPulseSequence(tmpStruct) ; 
            elseif isa(this,'mig.SIRqMTPulseSequence') 
            else        % basic PulseSequence
                varargout{1} = mig.PulseSequence(tmpStruct) ; 
            end
        elseif length(s) == 2 && strcmp(s(2).type,'.')
            varargout{1} = this.(s(2).subs)(s(1).subs{1});
        elseif length(s) == 3 && strcmp(s(2).type,'.') && strcmp(s(3).type,'()')
            % Implement this(indices).PropertyName(indices)
            ...
        else
        % Use built-in for any other expression
        [varargout{1:nargout}] = builtin('subsref',this,s);
        end
    case '{}'
        if length(s) == 1
%             varargout{1} = eval(class(this));            
%             mc = metaclass(varargout{1});
%             is_Constant = [mc.PropertyList.Constant];
%             is_Dependent = [mc.PropertyList.Dependent];
%             all_fields = {mc.PropertyList.Name};
%             fields = all_fields(~is_Constant & ~is_Dependent);
%             for i = 1:length(fields)
%                 varargout{1}.(fields{i}) = this.(fields{i})(s.subs{:},:);
%             end
        elseif length(s) == 2 && strcmp(s(2).type,'.')
            % Implement this{indices}.PropertyName
            ...
        else
        % Use built-in for any other expression
        [varargout{1:nargout}] = builtin('subsref',this,s);
        end
    otherwise
        error('Not a valid indexing expression')
end