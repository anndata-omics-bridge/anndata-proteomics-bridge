-- Rewrite local Markdown links (foo.md, foo.md#anchor) to .html so the
-- generated HTML navigates between rendered pages. Absolute URLs (with a
-- scheme like http:) and pure anchors are left untouched. The .md sources keep
-- their .md links, which is correct for GitHub / VS Code rendering.
function Link(el)
  if not el.target:match("^%a[%w+.-]*:") then
    el.target = el.target:gsub("%.md$", ".html"):gsub("%.md#", ".html#")
  end
  return el
end
