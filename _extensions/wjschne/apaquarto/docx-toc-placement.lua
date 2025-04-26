-- docx-toc-placement.lua
function Div(div)
  if div.identifier == "toc" and FORMAT:match 'docx' then
    -- Create a manual TOC field at this specific location
    local toc_field = pandoc.RawBlock('openxml', 
      '<w:p><w:pPr><w:pStyle w:val="Heading5"/><w:jc w:val="center"/></w:pPr>' ..
      '<w:r><w:t>Table of Contents</w:t></w:r></w:p>' ..
      '<w:p><w:r><w:fldChar w:fldCharType="begin"/></w:r>' ..
      '<w:r><w:instrText xml:space="preserve"> TOC \\o "1-3" \\h \\z \\u </w:instrText></w:r>' ..
      '<w:r><w:fldChar w:fldCharType="separate"/></w:r>' ..
      '<w:r><w:t>Contents will appear here when updated in Word</w:t></w:r>' ..
      '<w:r><w:fldChar w:fldCharType="end"/></w:r></w:p>')
    
    return toc_field
  end
  return div
end

-- Prevent automatic TOC generation
function Meta(meta)
  if FORMAT:match 'docx' then
    -- Ensure header-includes exists
    if not meta['header-includes'] then
      meta['header-includes'] = pandoc.MetaList{}
    end
    
    -- This empty comment helps override Quarto's automatic TOC insertion
    table.insert(
      meta['header-includes'],
      pandoc.MetaBlocks{
        pandoc.RawBlock('openxml', '<!-- No automatic TOC -->')
      }
    )
  end
  return meta
end

return {
  {
    Meta = Meta,
    Div = Div
  }
}