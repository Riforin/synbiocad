from django import forms
from django.utils.translation import ugettext as _
from django.utils.safestring import mark_safe
from django.forms.widgets import RadioSelect

class IndexForm(forms.Form):
	CHOICES = [['1','Insert DNA into a characterized locus'],
		       ['2','Edit an existing gene, e.g., delete or replace'],
		  	   ['3','Build a custom cassette']]

	choices = forms.ChoiceField(required=True, 
								widget = forms.RadioSelect, 
								choices = CHOICES,
								label = False
								)

class CharLocusForm(forms.Form):
	cutsite = forms.CharField(max_length = 10)

	CHOICES = [['1', 'Insert a pre-built custom cassette'],
			   ['2', 'Construct a cassette using standard promoters']]

	choices = forms.ChoiceField(required = True,
								widget = forms.RadioSelect,
								choices = CHOICES, 
								label = 'What would you like to do to this cutsite?')

class CharLocusPrebuiltForm(forms.Form):
	name = forms.CharField(max_length = 100,
						   label = 'What is the name of the ORF?:')
	sequence = forms.CharField(max_length = 1000000,
							   label = 'What is the sequence of the ORF?:')
