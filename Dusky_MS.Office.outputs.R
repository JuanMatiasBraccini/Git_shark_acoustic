library(ReporteRs)

#function for creating word tables
fn.word.table=function(WD,TBL,Doc.nm,caption,paragph,HdR.col,HdR.bg,Hdr.fnt.sze,Hdr.bld,
                       body.fnt.sze,Zebra,Zebra.col,Grid.col,Fnt.hdr,Fnt.body)
{
  mydoc = docx(Doc.nm)  #create r object
  
  # add title
  if(!is.na(caption))mydoc = addParagraph(mydoc, caption, stylename = "TitleDoc" )
  
  # add a paragraph
  if(!is.na(paragph))mydoc = addParagraph(mydoc , paragph, stylename="Citationintense")
  
  #add table
  MyFTable=FlexTable(TBL,add.rownames =F,
     header.cell.props = cellProperties(background.color=HdR.bg), 
     header.text.props = textProperties(color=HdR.col,font.size=Hdr.fnt.sze,
                        font.weight="bold",font.family =Fnt.hdr), 
     body.text.props = textProperties(font.size=body.fnt.sze,font.family =Fnt.body))
  
    # zebra stripes - alternate colored backgrounds on table rows
  if(Zebra=="YES") MyFTable = setZebraStyle(MyFTable, odd = Zebra.col, even = "white" )
  
  # table borders
  MyFTable = setFlexTableBorders(MyFTable,
      inner.vertical = borderNone(),inner.horizontal = borderNone(),
      outer.vertical = borderNone(),
      outer.horizontal = borderProperties(color=Grid.col, style="solid", width=4))
  
  # set columns widths (in inches)
  #MyFTable = setFlexTableWidths( MyFTable, widths = Col.width)
  
  
  
  mydoc = addFlexTable( mydoc, MyFTable) 
  # write the doc 
  writeDoc( mydoc, file = paste(Doc.nm,".docx",sep=''))
}


#function for creating word figures
fn.word.figure=function(WD,fig,Doc.nm,caption,paragph)
{
  options( "ReporteRs-fontsize" = 12 )  #set font
  mydoc = docx(Doc.nm)  #create r object
  
  # add title
  doc = addParagraph(mydoc, caption, stylename = "TitleDoc" )
  
  # add a paragraph
  doc = addParagraph(mydoc , paragph, stylename="Citationintense")
  
  # add a plot into mydoc 
  mydoc = addPlot( mydoc, fig)
  
  # write the doc 
  writeDoc( mydoc, file = paste(Doc.nm,".docx",sep=''))
}
