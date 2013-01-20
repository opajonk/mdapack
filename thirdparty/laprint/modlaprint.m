function modlaprint(figno,dirname,filename,varargin)
% MODLAPRINT Modifying the laprint-produced file for easier inclusion. See the
%   laprint documentation for more information. This script is run
%   in a similar way as laprint itself and run instead of laprint.
%
% Usage:
%   >>modlaprint(<figno>,<dirname>,<filename>,<opt1>,<optVal1>,...)
%
%   where
%   <figno>   (mandatory) is the MATLAB figure number,
%   <dirname> (mandatory) is the main directory for EPS and TeX
%   <filename>(mandatory) is the laprint produced EPS and TeX files
%                         placed in dirname/images/
%   {<opt>,<optVal>} (optional) are the laprint supported
%   options as-well as the modlaprint supported ones. The options
%   are declared as {'option name','its corresponding value',...}.
%   Current supported options are every
%   laprint option and the options for this script:
%   'texoutput' with value 'texfile'     -- Optional filename where
%                                        -- the laprint EPS and
%                                        -- TeX files are included
%                                        -- in dirname/texfile
%                                        -- If file allready
%                                        -- exists, the text is
%                                        -- appended to the file.
%                                        -- If a path is given containing
%                                        -- directories, absolute path
%                                        -- names for the images will be
%                                        -- written in the TeX-file
%   'makeps' with value '1' or '0'       -- run latex on dirname/texfile
%                                        -- run dvips on dirname/texfile
%   'ylabelrot' with value angle deg     -- Rotates the ylabel with
%                                        -- -90 as uppward allignement
%   'section'   with value 'some string' -- Will inclue a section
%   'subsection with value 'some string' -- Will inclue a subsection
%
% Example:
%   >> modlaprint(1,'teximages','imgfile',caption','some text',
%      'texoutput','texprint','section','Results','subsection','Forces');
%
%   Produces laprint `images/image1.eps' and `images/image1.tex' files
%   of MATLAB figure(1), passes a caption string both to laprint and
%   modlaprint, produces the file `tex/imagefile.tex' with graphic
%   inclusion, section and subsection inclusions. Note that the
%   laprint produced files are here located in the directoy `images'
%   and the modlaprint produced file is located in `tex'. However,
%   'tex/imagefile.tex' will include relative references to the
%   directoty 'images' which will lead to a `filenotfound'-error.
%   This is simply handled by creating a symlink to ../images with
%   name 'images'. There are advantages in using relative paths, e.g.
%   more robust site independence.
%
%   To include the figure in the TeX document, which is the way this
%   script does when the <texoutput> argument is specified, do e.g.:
%
%    \begin{figure}[!htb]
%     \centering{\figtext     % \figtext has the fontsize in the
%                             % figure. It is custom to have the
%                             % fontsize a factor 0.8 of the body
%                             % text. Created in the preamble by:
%                             % \newcommand{\figtext}{footnotesize}
%
%     \input{images/filename} % This will include the TeX-file produced by
%                             % laprint containing the PSfrag commands
%     \includegraphics[width=0.5\textwidth]
%     {images/filname}        % This will include the EPS-file produced by
%                             % laprint
%     \textsf{\caption{Some caption}}
%    \end{figure}
%
%   Author: Fredrik Sahlin, fresah@ltu.se
%
%   See also LAPRINT.

%   $Version: 1.1 $  $Date: 2005/07/25 $

curr_dir    = pwd;
if ~isdir(dirname)
    disp(['The destination directory `',dirname,''' does not exist']);
else
    if ~isdir([dirname,'/','images'])
        cd(dirname);
        %cd(curr_dir);
        mkdir('images');
        cd(curr_dir);
    end

    lapropt     = {};
    section     = '';
    subsection  = '';
    captionStr  = '';
    texoutput   = '';
    imgpath     = 'images';
    ctrlapr     = 0;
    ylabelrot   = 0;
    absTexPath  = 0;

    %-------------------------------
    % Checking for options, sorting
    % out laprint and modlaprint
    % options
    %-------------------------------
    if ~isempty(varargin)
        for jj = 1:2:length(varargin)
            if     strcmp(varargin{jj},'texoutput');
                texoutput   = varargin{jj+1};
                % If texoutput contains '/' we use absolute paths,
                % otherwise we use relative paths in tex-file.
                testtex = texoutput;
                [s,f,t] = regexp(testtex,'^(.*)/[^/]+$');
                if s
                    absimgpath = testtex(t{1}(1):t{1}(2));
                    absTexPath = 1;
                end
            elseif strcmp(varargin{jj},'section');      section     = varargin{jj+1};
            elseif strcmp(varargin{jj},'subsection');   subsection  = varargin{jj+1};
            elseif strcmp(varargin{jj},'makeps');       makeps      = varargin{jj+1};
            elseif strcmp(varargin{jj},'ylabelrot');
                ylabelangle = varargin{jj+1};
                ylabelrot   = 1;
            elseif strcmp(varargin{jj},'caption');      captionStr  = varargin{jj+1};
            else % Not a known option; pass it to laprint
                ctrlapr          = ctrlapr + 1;
                lapropt{ctrlapr} = varargin{jj};
                ctrlapr          = ctrlapr + 1;
                lapropt{ctrlapr} = varargin{jj+1};
            end
        end
    end

    pstest = 1;
    if ~(exist('makeps'))
        pstest = 0;
    else
        pstest = makeps;
    end

    set(0,'defaulttextinterpreter','none');
    warning off MATLAB:conversionToLogical;
    cd([dirname,'/images']);
    %cd([curr_dir,'/images']);
    %if regexp(dirname,'^/.*$')
    %  cd([dirname,'/images']);
    %else
    %  cd([curr_dir,'/',dirname,'/images']);
    %end
    laprint(figno,filename,lapropt{1:end});
    fixPSlinestyle([filename,'.eps']);

    disp(['Writing image: images/',filename])
    
    copyfile([filename,'.tex'],[filename,'.tmp']);
    MLTEXR = fopen([filename,'.tmp'],'r');
    MLTEXW = fopen([,filename,'.tex'],'w');
    cd(curr_dir);
    while 1
        tline = fgets(MLTEXR);
        if ~ischar(tline), break, end
        tline = regexprep(tline, '^(\\begin{psfrags}.*)$', ['% Modified by MODLAPRINT',char(10),'%>> $1']);
        tline = regexprep(tline, '^(\\parbox.*)$',         '%>> $1');
        tline = regexprep(tline, '^(\\resizeb.*)$',        '%>> $1');
        tline = regexprep(tline, '^(\\caption.*)$',        '%>> $1');
        tline = regexprep(tline, '^(\\label{fig.*)$',      '%>> $1');
        tline = regexprep(tline, '^(}%.*)$',               '%>> $1');
        tline = regexprep(tline, '^(\\end{psfrags}.*)$',   '%>> $1');
        tline = regexprep(tline, '^(\\psfrag{s01})\[t\]\[t\](.*)$', ...
            '$1[t][b][1]$2');
        if ylabelrot
            tline = regexprep(tline, '^(\\psfrag{s02})\[b\]\[b\](.*)$', ...
                ['$1[tr][tr][1][',ylabelangle,']$2']);
        end
        fprintf(MLTEXW,'%s',tline);
    end
    fclose(MLTEXR);
    fclose(MLTEXW);
    delete([dirname,'/images/',filename,'.tmp']);

    if texoutput
        texpath = [dirname,'/',texoutput];
        if absTexPath
            imgpath = [curr_dir,'/',dirname,'/images'];
            texpath = texoutput;
        end
        [stat,mess,messid]=fileattrib([texpath, '.tex']);
        if stat
            % If file exists, remove \end{document}
            copyfile([texpath,'.tex'],[texpath,'.tmp']);
            MLTEXR = fopen([texpath, '.tmp'],'r');
            MLTEXW = fopen([texpath, '.tex'],'w');
            while 1
                tline = fgets(MLTEXR);
                if ~ischar(tline), break, end
                tline = regexprep(tline, '^(\\end{document}.*)$', '');
                fprintf(MLTEXW,'%s',tline);
            end
            fclose(MLTEXR);
            fclose(MLTEXW);
            delete([texpath,'.tmp']);
        else
            % If file doesn't exists, write preamble
            MLTEXW = fopen([texpath '.tex'],'w');
            fprintf(MLTEXW,'%s\n','\documentclass[a4paper, twocolumn]{article}');
            fprintf(MLTEXW,'%s\n','\usepackage[latin1]{inputenc}');
            fprintf(MLTEXW,'%s\n','\usepackage[english]{babel}  ');
            fprintf(MLTEXW,'%s\n','\usepackage[dvips]{graphicx,color,psfrag} ');
            fprintf(MLTEXW,'%s\n','\usepackage{amssymb, amsmath}');
            fprintf(MLTEXW,'%s\n','\newcommand{\figtext}{\footnotesize}   % The font size for figures');
            fprintf(MLTEXW,'%s\n','\begin{document}');
            fprintf(MLTEXW,'%s\n','\selectlanguage{english}');
            fclose(MLTEXW);
        end
        %---------------------------------
        % Writing our own TeX-file and
        % including the laprint TeX and eps
        %---------------------------------

        % First, check if the TeX-file allready contains
        % the image filename to be apended. If true,
        % do not append
        MLTEX    = fopen([texpath '.tex'],'r');
        isImg = 0;
        while 1
            tline = fgetl(MLTEX);
            if tline == -1; break; end;
            tline = regexprep(tline,'^.*{(.*\.eps)}.*$','$1');
            if strcmp(tline, [imgpath,'/',filename,'.eps']);
                isImg = 1;
                break;
            end
        end
        fclose(MLTEX);
        if isImg
            disp('Image name exists in file,');
            disp(['Skipping appending to:  ',texpath, '.tex'])
        else
            disp(['Appending to:  ',texpath, '.tex'])
            MLTEX = fopen([texpath '.tex'],'a');
            if section
                fprintf(MLTEX,'%s\n','\cleardoublepage');
                fprintf(MLTEX,'%s\n\n',['\section{' section '}']);
            end
            if subsection
                fprintf(MLTEX,'%s\n\n',['\subsection{' subsection '}']);
            end
            fprintf(MLTEX,'%s\n', '\clearpage');
            fprintf(MLTEX,'%s\n', '\begin{figure}[!h]');
            fprintf(MLTEX,'%s\n', '  \centering{\figtext');
            fprintf(MLTEX,'%s\n',['  \input{' [imgpath,'/',filename,'.tex'] '}']);
            fprintf(MLTEX,'%s\n', '  \includegraphics[width=\textwidth]');
            fprintf(MLTEX,'%s\n',['  {' [imgpath,'/',filename,'.eps'] '}}']);
            if captionStr
                fprintf(MLTEX,'%s\n',['  \caption{',captionStr,'}']);
                fprintf(MLTEX,'%s\n',['  \label{fig:',filename,'}']);
            end
            fprintf(MLTEX,'%s\n\n','\end{figure}');
            fclose(MLTEX);
        end
        MLTEX = fopen([texpath '.tex'],'a');
        fprintf(MLTEX,'%s\n','\end{document}');
        fclose(MLTEX);
    end

    % Create a ps file if the option <makeps> is 1
    if pstest
        cd(dirname);
        system(['latex ',texoutput]);
        system(['dvips ',texoutput]);
        file    = dir([texoutput,'.*']);
        for ii = 1:length(file)
            dot = strfind(file(ii).name,'.');
            if ~strcmp(file(ii).name(dot:end),'.ps') & ~strcmp(file(ii).name(dot:end),'.tex')
                delete(file(ii).name);
            end
        end
        cd(curr_dir);
    end

    set(0,'defaulttextinterpreter','tex');
end