The structure of the wiki is more complex than need be which is due to github.

Since _Sidebar.md, _Footer.md and _Header.md are interpreted globally, the individual chapters of the wiki have gotten their own folders with their own _Sidebar.md and _Footer.md.

If you change something that affects the structure of a wiki page, i.e. adding a section, make sure to update the _Sidebar.md that belongs to that page.

## Linking content:

- wiki style: [[Text that is shown | Name-of-the-wiki-page]] which would link to Name-of-the-wiki-page.md anywhere within the github wiki repository

- markdown style: [Text that is shown](Name-of-the-wiki-page) or [Text that is shown](http://mylink.com)


## When creating a new Page.md:

- give a _meaningful_ and unique name! Use hyphens instead of underscores, e.g. deepTools-help-with-something.md as hyphens will not be shown in the /Pages section on the github wiki, thus making this index of all wiki pages more readable
- make sure to add the path to the file in the header of the text, e.g. [**HOME**](Home) > [**FAQ**](FAQs) > [New Page](My-new-wiki-page)
- copy the table of content to _Sidebar.md
