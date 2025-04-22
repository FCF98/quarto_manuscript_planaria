-- This filter prints the apa-note, if present
-- Do nothing if latex
if FORMAT == "latex" then
  return
end
-- Default word for note
local beginapanote = "Note"
-- Replace note word, if specified
local function getnote(m)
    if m.language and m.language["figure-table-note"] then
        beginapanote = pandoc.utils.stringify(m.language["figure-table-note"])
    end
end
local function apanote(elem)
  if elem.attributes["apa-note"] then
      hasnote = true
      -- If div contains another div with apa-note, do nothing
      elem.content:walk {
        Div = function(div) 
          if div.attributes["apa-note"] then
            hasnote = false
          end
        end
      }
      if hasnote then
        -- Make note with "Note" in bold but content not in italics
        local apanotepara = pandoc.Para({
          pandoc.Strong(pandoc.Str(beginapanote)), 
          pandoc.Str("."),
          pandoc.Space()
        })
        -- Add the note content without italics
        apanotepara.content:extend(quarto.utils.string_to_inlines(elem.attributes["apa-note"]))
        
        local apanote = pandoc.Div(apanotepara)
        
        apanote.attributes['custom-style'] = 'FigureNote'
        apanote.classes:extend({"FigureNote"})
        return {elem, apanote}
      end
    end
end
return {
  {Meta = getnote},
  {Div = apanote}
}