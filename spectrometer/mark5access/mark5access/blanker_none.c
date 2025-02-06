/***************************************************************************
 *   Copyright (C) 2007 by Walter Brisken                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
//===========================================================================
// SVN properties (DO NOT CHANGE)
//
// $Id: blanker_none.c,v 1.4 2012/02/21 13:35:44 jwagnerhki Exp $
// $HeadURL: https://svn.atnf.csiro.au/difx/libraries/mark5access/trunk/mark5access/blanker_none.c $
// $LastChangedRevision: 777 $
// $Author: jwagnerhki $
// $LastChangedDate: 2008-09-10 16:48:08 +0200 (Wed, 10 Sep 2008) $
//
//============================================================================

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "mark5access/mark5_stream.h"

int blanker_none(struct mark5_stream *ms)
{
	ms->log2blankzonesize = 30;
	ms->blankzonestartvalid[0] = 0;
	ms->blankzoneendvalid[0] = 1<<30;

	return 0;
}
