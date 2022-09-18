%%
%% This is file `mighelbib.bst',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% merlin.mbs  (with options: `head,exlang,lang,annote,nm-rev1,jnrlst,nmft,nmft-sc,nmand-rm,dt-end,yr-com,tit-qq,qt-g,qx,atit-u,vol-bf,volp-com,pgsep-s,jwdpg,pp-last,num-xser,numser,ser-vol,ser-ed,pre-edn,edpar,blk-com,spanish,pp,and-xcom,etal-xc,etal-it,revdata,eprint,url,url-blk,em-it,nfss,,{}')
%% ----------------------------------------
%% ***  ***
%% 
%% Copyright 1994-2005 Patrick W Daly
 % ===============================================================
 % IMPORTANT NOTICE:
 % This bibliographic style (bst) file has been generated from one or
 % more master bibliographic style (mbs) files, listed above.
 %
 % This generated file can be redistributed and/or modified under the terms
 % of the LaTeX Project Public License Distributed from CTAN
 % archives in directory macros/latex/base/lppl.txt; either
 % version 1 of the License, or any later version.
 % ===============================================================
 % Name and version information of the main mbs file:
 % \ProvidesFile{merlin.mbs}[2005/10/17 4.14 (PWD, AO, DPC)]
 %   For use with BibTeX version 0.99a or later
 %-------------------------------------------------------------------
 % This bibliography style file is intended for texts in
 %   SPANISH
 % This is a numerical citation style, and as such is standard LaTeX.
 % It requires no extra package to interface to the main text.
 % The form of the \bibitem entries is
 %   \bibitem{key}...
 % Usage of \cite is as follows:
 %   \cite{key} ==>>          [#]
 %   \cite[chap. 2]{key} ==>> [#, chap. 2]
 % where # is a number determined by the ordering in the reference list.
 % The order in the reference list is alphabetical by authors.
 %---------------------------------------------------------------------

ENTRY
  { address
    annote
    archive
    author
    booktitle
    chapter
    collaboration
    edition
    editor
    eid
    eprint
    howpublished
    institution
    journal
    key
    language
    month
    note
    number
    numpages
    organization
    pages
    publisher
    school
    series
    title
    type
    url
    volume
    year
  }
  {}
  { label }
INTEGERS { output.state before.all mid.sentence after.sentence after.block }
FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}
STRINGS { s t}
FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
        { add.period$ write$
          newline$
          "\newblock " write$
        }
        { output.state before.all =
            'write$
            { add.period$ " " * write$ }
          if$
        }
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}
FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}
FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}
FUNCTION {fin.entry}
{ add.period$
  write$
  newline$
  annote missing$
    { "\bibAnnoteFile{" cite$ * "}" * write$ newline$ }
    { "\bibAnnote{" cite$ * "}{" * annote * "}" * write$ newline$
    }
  if$
}
FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}
FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
        'skip$
        { after.sentence 'output.state := }
      if$
    }
  if$
}
FUNCTION {add.blank}
{  " " * before.all 'output.state :=
}
FUNCTION {date.block}
{
  skip$
}
FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}
FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}
FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}
FUNCTION {non.stop}
{ duplicate$
   "}" * add.period$
   #-1 #1 substring$ "." =
}
FUNCTION {new.block.checka}
{ empty$
    'skip$
    'new.block
  if$
}
FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}
FUNCTION {new.sentence.checka}
{ empty$
    'skip$
    'new.sentence
  if$
}
FUNCTION {new.sentence.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.sentence
  if$
}
FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}
FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "\textit{" swap$ * "}" * }
  if$
}
FUNCTION {bolden}
{ duplicate$ empty$
    { pop$ "" }
    { "\textbf{" swap$ * "}" * }
  if$
}
FUNCTION {scaps}
{ duplicate$ empty$
    { pop$ "" }
    { "\textsc{" swap$ * "}" * }
  if$
}
FUNCTION {bib.name.font}
{ scaps }
FUNCTION {bib.fname.font}
{ bib.name.font }
FUNCTION {tie.or.space.prefix}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$
}
FUNCTION {capitalize}
{ "u" change.case$ "t" change.case$ }
FUNCTION {space.word}
{ " " swap$ * " " * }
 % Here are the language-specific definitions for explicit words.
 % Each function has a name bbl.xxx where xxx is the English word.
 %-------------------------------------------------------------------
 % Begin module:
 % \ProvidesFile{spanish.mbs}[2003/11/06 4.2 (PWD/TB)]
 % The language selected here is SPANISH
FUNCTION {bbl.and}
{ "y"}
FUNCTION {bbl.etal}
{ "et~al." }
FUNCTION {bbl.editors}
{ "editores" }
FUNCTION {bbl.editor}
{ "editor" }
FUNCTION {bbl.edby}
{ "editado por" }
FUNCTION {bbl.edition}
{ "edici\'on" }
FUNCTION {bbl.volume}
{ "tomo" }
FUNCTION {bbl.of}
{ "de" }
FUNCTION {bbl.number}
{ "n\'umero" }
FUNCTION {bbl.nr}
{ "n\raise1ex\hbox{\underbar{\scriptsize o}}" }
FUNCTION {bbl.in}
{ "en" }
FUNCTION {bbl.pages}
{ "p\'ags." }
FUNCTION {bbl.page}
{ "p\'ag." }
FUNCTION {bbl.eidpp}
{ "p\'aginas" }
FUNCTION {bbl.chapter}
{ "cap\'{\i}tulo" }
FUNCTION {bbl.techrep}
{ "Informe T\'ecnico" }
FUNCTION {bbl.mthesis}
{ "Tesis de Maestría" }
FUNCTION {bbl.phdthesis}
{ "Tesis Doctoral" }
FUNCTION {bbl.first}
{ "primera" }
FUNCTION {bbl.second}
{ "segunda" }
FUNCTION {bbl.third}
{ "tercera" }
FUNCTION {bbl.fourth}
{ "cuarta" }
FUNCTION {bbl.fifth}
{ "quinta" }
FUNCTION {bbl.th}
{ "\raise1ex\hbox{\underbar{\scriptsize a}}" }
MACRO {jan} {"enero"}
MACRO {feb} {"febrero"}
MACRO {mar} {"marzo"}
MACRO {apr} {"abril"}
MACRO {may} {"mayo"}
MACRO {jun} {"junio"}
MACRO {jul} {"julio"}
MACRO {aug} {"agosto"}
MACRO {sep} {"septiembre"}
MACRO {oct} {"octubre"}
MACRO {nov} {"noviembre"}
MACRO {dec} {"diciembre"}
% End module: spanish.mbs

FUNCTION {bibinfo.check}
{ swap$
  duplicate$ missing$
    {
      pop$ pop$
      ""
    }
    { duplicate$ empty$
        {
          swap$ pop$
        }
        { swap$
          pop$
        }
      if$
    }
  if$
}
FUNCTION {bibinfo.warn}
{ swap$
  duplicate$ missing$
    {
      swap$ "missing " swap$ * " in " * cite$ * warning$ pop$
      ""
    }
    { duplicate$ empty$
        {
          swap$ "empty " swap$ * " in " * cite$ * warning$
        }
        { swap$
          pop$
        }
      if$
    }
  if$
}
FUNCTION {format.eprint}
{ eprint duplicate$ empty$
    'skip$
    { "\eprint"
      archive empty$
        'skip$
        { "[" * archive * "]" * }
      if$
      "{" * swap$ * "}" *
    }
  if$
}
FUNCTION {format.url}
{ url empty$
    { "" }
    { "\urlprefix\url{" url * "}" * }
  if$
}

STRINGS  { bibinfo}
INTEGERS { nameptr namesleft numnames }

FUNCTION {format.names}
{ 'bibinfo :=
  duplicate$ empty$ 'skip$ {
  's :=
  "" 't :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      duplicate$ #1 >
        { "{f.~}{vv~}{ll}{, jj}" }
        { "{vv~}{ll}{, f.}{, jj}" }
      if$
      format.name$
      bib.name.font
      bibinfo bibinfo.check
      't :=
      nameptr #1 >
        {
          namesleft #1 >
            { ", " * t * }
            {
              s nameptr "{ll}" format.name$ duplicate$ "others" =
                { 't := }
                { pop$ }
              if$
              t "others" =
                {
                  " " * bbl.etal emphasize *
                }
                {
                  bbl.and
                  space.word * t *
                }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
  } if$
}
FUNCTION {format.names.ed}
{
  'bibinfo :=
  duplicate$ empty$ 'skip$ {
  's :=
  "" 't :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      "{f.~}{vv~}{ll}{, jj}"
      format.name$
      bibinfo bibinfo.check
      't :=
      nameptr #1 >
        {
          namesleft #1 >
            { ", " * t * }
            {
              s nameptr "{ll}" format.name$ duplicate$ "others" =
                { 't := }
                { pop$ }
              if$
              t "others" =
                {
                  " " * bbl.etal emphasize *
                }
                {
                  bbl.and
                  space.word * t *
                }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
  } if$
}
FUNCTION {format.authors}
{ author "author" format.names
    duplicate$ empty$ 'skip$
    { collaboration "collaboration" bibinfo.check
      duplicate$ empty$ 'skip$
        { " (" swap$ * ")" * }
      if$
      *
    }
  if$
}
FUNCTION {get.bbl.editor}
{ editor num.names$ #1 > 'bbl.editors 'bbl.editor if$ }
FUNCTION {format.editors}
{ editor "editor" format.names duplicate$ empty$ 'skip$
    {
      " " *
      get.bbl.editor
   "(" swap$ * ")" *
      *
    }
  if$
}
FUNCTION {select.language}
{ duplicate$ empty$
    'pop$
    { language empty$
        'skip$
        { "{\selectlanguage{" language * "}" * swap$ * "}" * }
      if$
    }
    if$
}
FUNCTION {format.note}
{
 note empty$
    { "" }
    { note #1 #1 substring$
      duplicate$ "{" =
        'skip$
        { output.state mid.sentence =
          { "l" }
          { "u" }
        if$
        change.case$
        }
      if$
      note #2 global.max$ substring$ * "note" bibinfo.check
    }
  if$
}
FUNCTION {format.title}
{ title
  "title" bibinfo.check
  duplicate$ empty$ 'skip$
    {
      "\enquote{" swap$ *
      "}, " *
      select.language
    }
  if$
}
FUNCTION {end.quote.title}
{ title empty$
    'skip$
    { before.all 'output.state := }
  if$
}
FUNCTION {output.bibitem}
{ newline$
  "\bibitem{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}
FUNCTION {if.digit}
{ duplicate$ "0" =
  swap$ duplicate$ "1" =
  swap$ duplicate$ "2" =
  swap$ duplicate$ "3" =
  swap$ duplicate$ "4" =
  swap$ duplicate$ "5" =
  swap$ duplicate$ "6" =
  swap$ duplicate$ "7" =
  swap$ duplicate$ "8" =
  swap$ "9" = or or or or or or or or or
}
FUNCTION {n.separate}
{ 't :=
  ""
  #0 'numnames :=
  { t empty$ not }
  { t #-1 #1 substring$ if.digit
      { numnames #1 + 'numnames := }
      { #0 'numnames := }
    if$
    t #-1 #1 substring$ swap$ *
    t #-2 global.max$ substring$ 't :=
    numnames #5 =
      { duplicate$ #1 #2 substring$ swap$
        #3 global.max$ substring$
        "\," swap$ * *
      }
      'skip$
    if$
  }
  while$
}
FUNCTION {n.dashify}
{
  n.separate
  't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
        { t #1 #2 substring$ "--" = not
            { "--" *
              t #2 global.max$ substring$ 't :=
            }
            {   { t #1 #1 substring$ "-" = }
                { "-" *
                  t #2 global.max$ substring$ 't :=
                }
              while$
            }
          if$
        }
        { t #1 #1 substring$ *
          t #2 global.max$ substring$ 't :=
        }
      if$
    }
  while$
}
FUNCTION {word.in}
{ bbl.in
  " " * }
FUNCTION {format.date}
{
  month "month" bibinfo.check
  duplicate$ empty$
  year  "year"  bibinfo.check duplicate$ empty$
    { swap$ 'skip$
        { "there's a month but no year in " cite$ * warning$ }
      if$
      *
    }
    { swap$ 'skip$
        {
          swap$
          " " * swap$
        }
      if$
      *
    }
  if$
  duplicate$ empty$
    'skip$
    {
      before.all 'output.state :=
    ", " swap$ *
    }
  if$
}
FUNCTION {format.btitle}
{ title "title" bibinfo.check
  duplicate$ empty$ 'skip$
    {
      emphasize
      select.language
    }
  if$
}
FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}
FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { bbl.volume volume tie.or.space.prefix
      "volume" bibinfo.check * *
      series "series" bibinfo.check
      duplicate$ empty$ 'pop$
        { emphasize ", " * swap$ * }
      if$
      "volume and number" number either.or.check
    }
  if$
}
FUNCTION {format.number.series}
{ volume empty$
    { number empty$
        { series field.or.null }
        { series empty$
            { number "number" bibinfo.check }
            { output.state mid.sentence =
                { bbl.number }
                { bbl.number capitalize }
              if$
              number tie.or.space.prefix "number" bibinfo.check * *
              bbl.in space.word *
              series "series" bibinfo.check *
            }
          if$
        }
      if$
    }
    { "" }
  if$
}
FUNCTION {is.num}
{ chr.to.int$
  duplicate$ "0" chr.to.int$ < not
  swap$ "9" chr.to.int$ > not and
}
FUNCTION {extract.num}
{ duplicate$ 't :=
  "" 's :=
  { t empty$ not }
  { t #1 #1 substring$
    t #2 global.max$ substring$ 't :=
    duplicate$ is.num
      { s swap$ * 's := }
      { pop$ "" 't := }
    if$
  }
  while$
  s empty$
    'skip$
    { pop$ s }
  if$
}
FUNCTION {convert.edition}
{ extract.num "l" change.case$ 's :=
  s "first" = s "1" = or
    { bbl.first 't := }
    { s "second" = s "2" = or
        { bbl.second 't := }
        { s "third" = s "3" = or
            { bbl.third 't := }
            { s "fourth" = s "4" = or
                { bbl.fourth 't := }
                { s "fifth" = s "5" = or
                    { bbl.fifth 't := }
                    { s #1 #1 substring$ is.num
                        { s bbl.th * 't := }
                        { edition 't := }
                      if$
                    }
                  if$
                }
              if$
            }
          if$
        }
      if$
    }
  if$
  t
}
FUNCTION {format.edition}
{ edition duplicate$ empty$ 'skip$
    {
      convert.edition
      output.state mid.sentence =
        { "l" }
        { "t" }
      if$ change.case$
      "edition" bibinfo.check
      " " * bbl.edition *
    }
  if$
}
INTEGERS { multiresult }
FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
        { #1 'multiresult := }
        { t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}
FUNCTION {format.pages}
{ pages duplicate$ empty$ 'skip$
    { duplicate$ multi.page.check
        {
          bbl.pages swap$
          n.dashify
        }
        {
          bbl.page swap$
        }
      if$
      tie.or.space.prefix
      "pages" bibinfo.check
      * *
    }
  if$
}
FUNCTION {format.journal.pages}
{ pages duplicate$ empty$ 'pop$
    { swap$ duplicate$ empty$
        { pop$ pop$ format.pages }
        {
          ", " *
          swap$
          n.dashify
          pages multi.page.check
            'bbl.pages
            'bbl.page
          if$
          swap$ tie.or.space.prefix
          "pages" bibinfo.check
          * *
          *
        }
      if$
    }
  if$
}
FUNCTION {format.journal.eid}
{ eid "eid" bibinfo.check
  duplicate$ empty$ 'pop$
    { swap$ duplicate$ empty$ 'skip$
      {
          ", " *
      }
      if$
      swap$ *
      numpages empty$ 'skip$
        { bbl.eidpp numpages tie.or.space.prefix
          "numpages" bibinfo.check * *
          " (" swap$ * ")" * *
        }
      if$
    }
  if$
}
FUNCTION {format.vol.num.pages}
{ volume field.or.null
  duplicate$ empty$ 'skip$
    {
      "volume" bibinfo.check
    }
  if$
  bolden
  number "number" bibinfo.check duplicate$ empty$ 'skip$
    {
      swap$ duplicate$ empty$
        { "there's a number but no volume in " cite$ * warning$ }
        'skip$
      if$
      swap$
      "(" swap$ * ")" *
    }
  if$ *
}
FUNCTION {format.chapter.pages}
{ chapter empty$
    { "" }
    { type empty$
        { bbl.chapter }
        { type "l" change.case$
          "type" bibinfo.check
        }
      if$
      chapter tie.or.space.prefix
      "chapter" bibinfo.check
      * *
    }
  if$
}
FUNCTION {format.booktitle}
{
  booktitle "booktitle" bibinfo.check
  emphasize
}
FUNCTION {format.in.ed.booktitle}
{ format.booktitle duplicate$ empty$ 'skip$
    {
      format.bvolume duplicate$ empty$ 'pop$
        { ", " swap$ * * }
      if$
      editor "editor" format.names.ed duplicate$ empty$ 'pop$
        {
          " " *
          get.bbl.editor
          "(" swap$ * "), " *
          * swap$
          * }
      if$
      word.in swap$ *
    }
  if$
}
FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  month empty$ year empty$ note empty$
  and and and and and
  key empty$ not and
    { "all relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}
FUNCTION {format.thesis.type}
{ type duplicate$ empty$
    'pop$
    { swap$ pop$
      "t" change.case$ "type" bibinfo.check
    }
  if$
}
FUNCTION {format.tr.number}
{ number "number" bibinfo.check
  type duplicate$ empty$
    { pop$ bbl.techrep }
    'skip$
  if$
  "type" bibinfo.check
  swap$ duplicate$ empty$
    { pop$ "t" change.case$ }
    { tie.or.space.prefix * * }
  if$
}
FUNCTION {format.article.crossref}
{
  key duplicate$ empty$
    { pop$
      journal duplicate$ empty$
        { "need key or journal for " cite$ * " to crossref " * crossref * warning$ }
        { "journal" bibinfo.check emphasize word.in swap$ * }
      if$
    }
    { word.in swap$ * " " *}
  if$
  " \cite{" * crossref * "}" *
}
FUNCTION {format.crossref.editor}
{ editor #1 "{vv~}{ll}" format.name$
  "editor" bibinfo.check
  editor num.names$ duplicate$
  #2 >
    { pop$
      "editor" bibinfo.check
      " " * bbl.etal
      emphasize
      *
    }
    { #2 <
        'skip$
        { editor #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
            {
              "editor" bibinfo.check
              " " * bbl.etal
              emphasize
              *
            }
            {
             bbl.and space.word
              * editor #2 "{vv~}{ll}" format.name$
              "editor" bibinfo.check
              *
            }
          if$
        }
      if$
    }
  if$
}
FUNCTION {format.book.crossref}
{ volume duplicate$ empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      pop$ word.in
    }
    { bbl.volume
      swap$ tie.or.space.prefix "volume" bibinfo.check * * bbl.of space.word *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
        { series empty$
            { "need editor, key, or series for " cite$ * " to crossref " *
              crossref * warning$
              "" *
            }
            { series emphasize * }
          if$
        }
        { key * }
      if$
    }
    { format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}
FUNCTION {format.incoll.inproc.crossref}
{
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
        { format.booktitle duplicate$ empty$
            { "need editor, key, or booktitle for " cite$ * " to crossref " *
              crossref * warning$
            }
            { word.in swap$ * }
          if$
        }
        { word.in key * " " *}
      if$
    }
    { word.in format.crossref.editor * " " *}
  if$
  " \cite{" * crossref * "}" *
}
FUNCTION {format.org.or.pub}
{ 't :=
  ""
  address empty$ t empty$ and
    'skip$
    {
      t empty$
        { address "address" bibinfo.check *
        }
        { t *
          address empty$
            'skip$
            { ", " * address "address" bibinfo.check * }
          if$
        }
      if$
    }
  if$
}
FUNCTION {format.publisher.address}
{ publisher "publisher" bibinfo.warn format.org.or.pub
}
FUNCTION {format.organization.address}
{ organization "organization" bibinfo.check format.org.or.pub
}
FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check
  end.quote.title
  crossref missing$
    {
      journal
      "journal" bibinfo.check
      emphasize
      "journal" output.check
      format.vol.num.pages output
    }
    { format.article.crossref output.nonnull
    }
  if$
  eid empty$
    { format.journal.pages }
    { format.journal.eid }
  if$
  format.note output
  format.date "year" output.check
  format.eprint output
  format.url output
  fin.entry
}
FUNCTION {book}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
    }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  format.btitle "title" output.check
  crossref missing$
    { format.bvolume output
      format.edition output
      format.number.series output
      format.publisher.address output
    }
    {
      format.book.crossref output.nonnull
    }
  if$
  format.note output
  format.date "year" output.check
  format.eprint output
  format.url output
  fin.entry
}
FUNCTION {booklet}
{ output.bibitem
  format.authors output
  format.title "title" output.check
  end.quote.title
  howpublished "howpublished" bibinfo.check output
  address "address" bibinfo.check output
  format.note output
  format.date output
  format.eprint output
  format.url output
  fin.entry
}
FUNCTION {inbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
    }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  format.btitle "title" output.check
  crossref missing$
    {
      format.bvolume output
      format.chapter.pages "chapter and pages" output.check
      format.edition output
      format.number.series output
      format.publisher.address output
    }
    {
      format.chapter.pages "chapter and pages" output.check
      format.book.crossref output.nonnull
    }
  if$
  format.pages "pages" output.check
  format.note output
  format.date "year" output.check
  format.eprint output
  format.url output
  fin.entry
}
FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check
  end.quote.title
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.edition output
      format.chapter.pages output
      format.number.series output
      format.publisher.address output
    }
    { format.incoll.inproc.crossref output.nonnull
      format.chapter.pages output
    }
  if$
  format.pages "pages" output.check
  format.note output
  format.date "year" output.check
  format.eprint output
  format.url output
  fin.entry
}
FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check
  end.quote.title
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.number.series output
      publisher empty$
        { format.organization.address output }
        { organization "organization" bibinfo.check output
          format.publisher.address output
        }
      if$
    }
    { format.incoll.inproc.crossref output.nonnull
    }
  if$
  format.pages "pages" output.check
  format.note output
  format.date "year" output.check
  format.eprint output
  format.url output
  fin.entry
}
FUNCTION {conference} { inproceedings }
FUNCTION {manual}
{ output.bibitem
  author empty$
    { organization "organization" bibinfo.check
      duplicate$ empty$ 'pop$
        { output
          address "address" bibinfo.check output
        }
      if$
    }
    { format.authors output.nonnull }
  if$
  format.btitle "title" output.check
  author empty$
    { organization empty$
        {
          address "address" bibinfo.check output
        }
        'skip$
      if$
    }
    {
      organization "organization" bibinfo.check output
      address "address" bibinfo.check output
    }
  if$
  format.edition output
  format.note output
  format.date output
  format.eprint output
  format.url output
  fin.entry
}
FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  format.btitle
  "title" output.check
  bbl.mthesis format.thesis.type output.nonnull
  school "school" bibinfo.warn output
  address "address" bibinfo.check output
  format.note output
  format.date "year" output.check
  format.eprint output
  format.url output
  fin.entry
}
FUNCTION {misc}
{ output.bibitem
  format.authors output
  format.title output
  end.quote.title
  howpublished "howpublished" bibinfo.check output
  format.note output
  format.date output
  format.eprint output
  format.url output
  fin.entry
  empty.misc.check
}
FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  format.btitle
  "title" output.check
  bbl.phdthesis format.thesis.type output.nonnull
  school "school" bibinfo.warn output
  address "address" bibinfo.check output
  format.note output
  format.date "year" output.check
  format.eprint output
  format.url output
  fin.entry
}
FUNCTION {proceedings}
{ output.bibitem
  editor empty$
    { organization "organization" bibinfo.check output
    }
    { format.editors output.nonnull }
  if$
  format.btitle "title" output.check
  format.bvolume output
  editor empty$
    { publisher empty$
        'skip$
        {
          format.number.series output
          format.publisher.address output
        }
      if$
    }
    { publisher empty$
        {
          format.organization.address output }
        {
          organization "organization" bibinfo.check output
          format.publisher.address output
        }
      if$
     }
  if$
  format.note output
  format.date "year" output.check
  format.eprint output
  format.url output
  fin.entry
}
FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  format.title
  "title" output.check
  end.quote.title
  format.tr.number emphasize output.nonnull
  institution "institution" bibinfo.warn output
  address "address" bibinfo.check output
  format.note output
  format.date "year" output.check
  format.eprint output
  format.url output
  fin.entry
}
FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  format.title "title" output.check
  end.quote.title
  format.note "note" output.check
  format.date output
  format.eprint output
  format.url output
  fin.entry
}
FUNCTION {default.type} { misc }
READ
FUNCTION {sortify}
{ purify$
  "l" change.case$
}
INTEGERS { len }
FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}
FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      "{vv{ } }{ll{ }}{  f{ }}{  jj{ }}"
      format.name$ 't :=
      nameptr #1 >
        {
          "   "  *
          namesleft #1 = t "others" = and
            { "zzzzz" * }
            { t sortify * }
          if$
        }
        { t sortify * }
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}
FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}
FUNCTION {author.sort}
{ author empty$
    { key empty$
        { "to sort, need author or key in " cite$ * warning$
          ""
        }
        { key sortify }
      if$
    }
    { author sort.format.names }
  if$
}
FUNCTION {author.editor.sort}
{ author empty$
    { editor empty$
        { key empty$
            { "to sort, need author, editor, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { editor sort.format.names }
      if$
    }
    { author sort.format.names }
  if$
}
FUNCTION {author.organization.sort}
{ author empty$
    { organization empty$
        { key empty$
            { "to sort, need author, organization, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { "The " #4 organization chop.word sortify }
      if$
    }
    { author sort.format.names }
  if$
}
FUNCTION {editor.organization.sort}
{ editor empty$
    { organization empty$
        { key empty$
            { "to sort, need editor, organization, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { "The " #4 organization chop.word sortify }
      if$
    }
    { editor sort.format.names }
  if$
}
FUNCTION {presort}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.sort
    { type$ "proceedings" =
        'editor.organization.sort
        { type$ "manual" =
            'author.organization.sort
            'author.sort
          if$
        }
      if$
    }
  if$
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}
ITERATE {presort}
SORT
STRINGS { longest.label }
INTEGERS { number.label longest.label.width }
FUNCTION {initialize.longest.label}
{ "" 'longest.label :=
  #1 'number.label :=
  #0 'longest.label.width :=
}
FUNCTION {longest.label.pass}
{ number.label int.to.str$ 'label :=
  number.label #1 + 'number.label :=
  label width$ longest.label.width >
    { label 'longest.label :=
      label width$ 'longest.label.width :=
    }
    'skip$
  if$
}
EXECUTE {initialize.longest.label}
ITERATE {longest.label.pass}
FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{"  longest.label  * "}" *
  write$ newline$
  "\ProvideTextCommand{\guillemotleft}{OT1}{%"
  write$ newline$
  "  \leavevmode\raise .27ex\hbox{$\scriptscriptstyle\ll$}}"
  write$ newline$
  "\ProvideTextCommand{\guillemotright}{OT1}{%"
  write$ newline$
  "  \leavevmode\raise .27ex\hbox{$\scriptscriptstyle\gg$}}"
  write$ newline$
  "\newcommand{\enquote}[1]{\guillemotleft#1\guillemotright}"
  write$ newline$
  "\providecommand{\url}[1]{\texttt{#1}}"
  write$ newline$
  "\providecommand{\urlprefix}{URL }"
  write$ newline$
  "\providecommand{\selectlanguage}[1]{\relax}"
  write$ newline$
  "\providecommand{\bibAnnoteFile}[1]{%"
  write$ newline$
  "  \IfFileExists{#1}{\begin{quotation}\noindent\textsc{Key:} #1\\"
  write$ newline$
  "  \textsc{Annotation:}\ \input{#1}\end{quotation}}{}}"
  write$ newline$
  "\providecommand{\bibAnnote}[2]{%"
  write$ newline$
  "  \begin{quotation}\noindent\textsc{Key:} #1\\"
  write$ newline$
  "  \textsc{Annotation:}\ #2\end{quotation}}"
  write$ newline$
  "\providecommand{\eprint}[2][]{\url{#2}}"
  write$ newline$
}
EXECUTE {begin.bib}
EXECUTE {init.state.consts}
ITERATE {call.type$}
FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}
EXECUTE {end.bib}
%% End of customized bst file
%%
%% End of file `mighelbib.bst'.
