# ==============================================================================
.PHONY: all html latex clean clean-latex clean-html

html:    
	doxygen Doxyfile_html

dev:    
	doxygen Doxyfile_dev

latex:    
	doxygen Doxyfile_latex
	cd latex; make

all:    
	doxygen Doxyfile

clean:
	-rm -fr html/* latex/*

clean-html:
	-rm -fr latex/*

clean-latex:
	-rm -fr latex/*

# ==============================================================================
